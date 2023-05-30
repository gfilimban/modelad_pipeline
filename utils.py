import pysam
import pandas as pd
import re
import math
import gzip
import pyranges as pr
import numpy as np
import scanpy as sc
import swan_vis as swan

def process_meta(meta_fname):
    meta = pd.read_csv(meta_fname, sep='\t')
    meta['genotype'] = meta.genotype.str.replace(' ', '_')
    meta['genotype'] = meta.genotype.str.replace('/', '_')

    meta['age'] = meta.age.str.replace(' ', '_')

    dupe_ids = meta.loc[meta.mouse_id.duplicated(), 'mouse_id'].tolist()
    if len(dupe_ids) > 0:
        raise ValueError(f'Found duplicated mouse ids {dupe_ids} in mouse metadata')

    return meta

def parse_config_file(fname,
                      meta_fname,
                      datasets_per_run,
                      auto_dedupe=True):

    """
    Parameters:
        fname (str): Path to config file fname. One line per input fastq.
        meta_fname (str): Path to file with metadata information.
        datasets_per_run (int): Number of datasets to process in each TALON run
        auto_dedupe (bool): Automatically deduplicate duplicate fastqs that result from
            successive Porechop rounds

    Returns:
        df (pandas DataFrame): DF w/ pipeline information; one line per fastq
        dataset_df (pandas DataFrame): DF w/ dataset information; one line per mouse
    """

    df = pd.read_csv(fname, sep='\t')

    ############ Dataset + flowcell df

    # get flowcell
    exp = '.*\/[\w-]+_(\d+)(?:_t\d+)?\.fastq(?:.gz)?'
    df['flowcell'] = df.fname.str.extract(exp)

    # check to make sure the same file stem isn't there more than once
    # (can happen if different flow cells needed different amounts of chopping)
    # df['file_stem'] = df.basename.str.rsplit('_', n=1, expand=True)[0]
    exp = '.*\/([\w-]+_\d+)(?:_t\d+)?\.fastq(?:.gz)?'
    df['file_stem'] = df.fname.str.extract(exp)
    df['chop_num'] = df.basename.str.rsplit('.fastq', expand=True)[0].str.rsplit('_t', expand=True)[1].astype(float)
    if df.file_stem.duplicated().any():
        dupe_stems = df.loc[df.file_stem.duplicated(keep=False), 'basename'].tolist()
        if not auto_dedupe:
            raise ValueError(f'Files {dupe_stems} seem to be duplicated. Check config file.')
        else:
            print(f'Files {dupe_stems} seem to be duplicated. Automatically removing lower chop numbers')
            df = df.sort_values(by='chop_num', ascending=False)
            df = df.drop_duplicates(subset='file_stem', keep='first')

    # extract the sample name
    temp = df.basename.str.split('_', expand=True)[[0,1]]#.str.join('_')
    df['sample_temp'] = temp[0]+'_'+temp[1]

    # extract the mouse id
    df['mouse_id'] = df['sample_temp'].str.split('_', expand=True)[1]

    # extract the "study" name
    exp = '^(ad[0-9]+)'
    df['study'] = df.basename.str.extract(exp)


    # merge in metadata
    meta = process_meta(meta_fname)
    df['mouse_id'] = df['mouse_id'].astype('int')
    df = df.merge(meta, how='left', on='mouse_id')

    # get tech rep numbers -- each mouse has multiple reps
    # and are therefore technical reps
    df['flowcell'] = df.sort_values(['genotype', 'mouse_id'],
                                ascending=[True, True])\
                                .groupby(['mouse_id']) \
                                .cumcount() + 1

    # sample should be the genotype + age + sex + tissue
    df['sample'] = df.genotype+'_'+ \
                   df.sex+'_'+ \
                   df.age+'_'+ \
                   df.tissue

    # get biorep numbers -- each mouse_id is a different mouse
    # and therefore a different biorep
    temp = df[['sample', 'mouse_id']].drop_duplicates()
    temp.reset_index(inplace=True, drop=True)
    temp['biorep_num'] = temp.sort_values(['sample', 'mouse_id'],
                                ascending=[True, True])\
                                .groupby(['sample']) \
                                .cumcount()+1
    df = df.merge(temp, how='left',
                  on=['sample', 'mouse_id'])

    # talon dataset should be sample + bio rep
    df['dataset'] = df['sample']+'_'+df['biorep_num'].astype(str)

    # # dataset should be sample + bio rep + flow cel
    # df['dataset'] = df['talon_dataset']+'_'+df['flowcell'].astype(str)

    ############ TALON dataset df

    # create a dataset-level df that will represent the aggregate
    cols = ['sample', 'mouse_id', 'genotype', 'sex', 'study', \
            'age', 'tissue', 'biorep_num', 'dataset', 'platform']
    dataset_df = df[cols].drop_duplicates()

    # get the talon run number these will go into
    dataset_df = dataset_df.sort_values(by='study', ascending=False)
    dataset_df = dataset_df.reset_index(drop=True)
    dataset_df['talon_run_num'] = np.nan
    for ind, entry in dataset_df.iterrows():
        curr_study = entry.study
        # first entry
        if ind == 0:
            run_num = 0
            study_ind = 0
            prev_study = curr_study

        # when we find a new study
        if curr_study != prev_study:
            run_num = 0
            study_ind = 0

        # when we hit the max. # datasets / run
        if study_ind % datasets_per_run == 0:
            run_num += 1

        # actual assignment of number
        dataset_df.loc[ind, 'talon_run_num'] = run_num

        # keep track of last study that was assigned a number
        prev_study = curr_study

        # inc study ind
        study_ind += 1

    # clean up some typing
    dataset_df['talon_run_num'] = dataset_df.talon_run_num.astype(int)
    df['flowcell'] = df.flowcell.astype(str)

    return df, dataset_df

def rev_comp(seq):
    """ Returns the reverse complement of a DNA sequence,
        retaining the case of each letter"""
    complement = ""

    for base in seq:
        if base == "A": complement += "T"
        elif base == "T": complement += "A"
        elif base == "G": complement += "C"
        elif base == "C": complement += "G"
        elif base == "N": complement += "N"
        elif base == "a": complement += "t"
        elif base == "t": complement += "a"
        elif base == "g": complement += "c"
        elif base == "c": complement += "g"
        elif base == "n": complement += "n"
        elif base == "*": complement += "*"
        else:
            complement += base
            print("Warning: reverse complement function encountered unknown base " + "'" + base + "'")

    reverseComplement = complement[::-1]

    return reverseComplement

def flip_fastq(infile, outfile):

    reverse = ''

    if infile.endswith('.gz'):
        f = gzip.open(infile, 'rt')
    else:
        f = open(infile, 'r')
    o = open(outfile, 'w')

    for line in f:
        line = line.strip()
        # if line starts with @, write identifier to ouput file
        if re.match('^@.{19,49}', line):
            o.write(line + "\n")
        elif re.match('^[ATCG].*', line):
            seq = line
            # if there is a TTTTTT repeat at begining of read, reverse sequence then write to output file
            if re.match('TTTTTTTTTTTTTTTTTT*', seq):
                line = rev_comp(seq)
                o.write(line + '\n')
                reverse = True
            # if there is no TTTTT repeat at beginng of read, write to output file
            else:
                o.write(line + '\n')
                reverse = False
        else:
            if not re.match('^@ATCG([+]\n)', line):
                qual = line
                # if the sequence line was reversed, reverse the quality score line then write to ouput
                if reverse == True:
                    line = qual[::-1]
                    o.write(line + '\n')
                # if the sequence line was not reversed, write quality score line to output
                else:
                    o.write(line + '\n')
    o.close()
    f.close()

def reverse_alignment(infile, outfile, threads=1):

    reverse_strand = {0: 16, 16: 0}

    if infile.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    if outfile.endswith('bam'):
        out_mode = 'wb'
    else:
        out_mode = 'w'
    input =  pysam.AlignmentFile(infile, in_mode, threads=threads)
    output = pysam.AlignmentFile(outfile, out_mode, template=input, threads=threads)

    for read in input:
        if read.has_tag('ts') and read.flag in reverse_strand:
            if read.get_tag('ts') == '-':
                read.flag = reverse_strand[read.flag]
                read.set_tag('ts', '+')

        output.write(read)
    input.close()
    output.close()

################################################################################
############################ LAPA file post-proc ###############################
################################################################################

# def gtf_add_rescue_ism_cat(gtf):
#     """
#     Update LAPA GTF w/ ISM rescue category for those that were assigned new
#     ends from LAPA
#     """
#     gtf = pr.read_gtf(input.lapa_gtf, as_df=True)
#     gtf.loc[(gtf.transcript_id.str.contains('#'))&(gtf.ISM_transcript=='TRUE'),
#              'transcript_novelty'] = 'ISM_rescue'
#     gtf = pr.PyRanges(gtf)
#     return gtf

def ab_add_rescue_ism_cat(ab):
    """
    Update LAPA abundance w/ ISM rescue category for those that were assigned new
    ends from LAPA
    """
    df = pd.read_csv(ab, sep='\t')
    df.loc[(df.annot_transcript_id.str.contains('#'))&(df.transcript_novelty=='ISM'), 'transcript_novelty'] = 'ISM_rescue'
    return df

def filter_lapa_on_nov(df,
                        t_novs=['Known', 'NIC', 'NNC', 'ISM_rescue'],
                        g_novs=['Known']):
    """
    Filter LAPA output based on gene and transcript novelty.

    Input:
        df (pandas df): Abundance table from LAPA
        t_nov (list of str): Transcript novelty categories to keep
        g_nov (list of str): Gene novelty categories to keep

    Returns:
        filt_df (pandas df): DataFrame of gene id, transcript id passing filt
    """
    df = df.loc[df.transcript_novelty.isin(t_novs)]
    df = df.loc[df.gene_novelty.isin(g_novs)]
    filt_df = df[['annot_gene_id', 'annot_transcript_id']].drop_duplicates()
    filt_df = filt_df.rename({'annot_gene_id':'gid',
                              'annot_transcript_id': 'tid'}, axis=1)
    return filt_df

def filter_spikes(gtf):
    """
    Filter LAPA output based on SIRV / ERCC status

    Input:
        gtf (str): GTF path from LAPA
    Returns:
        filt_df (pandas df): DataFrame of gene id, transcript id passing filt
    """
    df = pr.read_gtf(input.lapa_gtf, as_df=True)
    df = df.loc[~df.Chromosome.str.contains('SIRV')]
    df = df.loc[~df.Chromosome.str.contains('ERCC')]
    df = df.locc[df.Feature == 'transcript']
    filt_df = df[['gene_id', 'transcript_id']].drop_duplicates()
    filt_df = filt_df.rename({'gene_id':'gid',
                              'transcript_id':'tid'}, axis=1)
    return filt_df

def get_ids_from_pass_list(filt_list):
    filt_df = pd.read_csv(filt_list, sep='\t')
    gids = filt_df.gid.unique().tolist()
    tids = filt_df.tid.unique().tolist()
    return gids, tids

def filt_lapa_ab(ab, filt_list):
    """
    Filter LAPA abundance using a TALON-style pass list
    """
    df = pd.read_csv(ab, sep='\t')
    gids, tids = get_ids_from_pass_list(filt_list)
    df = df.loc[(df.annot_gene_id.isin(gids))&(df.annot_transcript_id.isin(tids))]
    return df

def filt_lapa_gtf(gtf, filt_list):
    """
    Filter LAPA GTF using a TALON-style pass list
    """
    gtf = pr.read_gtf(gtf).as_df()
    gids, tids = get_ids_from_pass_list(filt_list)

    # first filter on tids
    gtf = gtf.loc[(gtf.transcript_id.isin(tids))|(gtf.Feature=='gene')]

    # then filter on gids
    gtf = gtf.loc[(gtf.gene_id.isin(gids))]

    gtf = pr.PyRanges(gtf)
    return gtf

def save_swan_adata(swan_file,
                    ofile,
                    how='iso'):
    """
    Save anndata obj from Swan.

    Input:
        swan_file (str): Input SwanGraph file
        ofile (str): Output AnnData file
        how (str): {'iso', 'tss', 'tes', 'ic', 'edge', 'loc', 'gene'}
    """
    sg = swan.read(swan_file)
    if how == 'gene':
        adata = sg.gene_adata
    elif how == 'iso':
        adata = sg.adata
    else:
        raise ValueError("You haven't implemented this yet.")
    adata.write(ofile)
