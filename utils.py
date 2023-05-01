import pysam
import pandas as pd
import re
import math
import gzip
import pyranges as pr
import numpy as np

def parse_config_file(fname,
                      datasets_per_run,
                      auto_dedupe=True):

    df = pd.read_csv(fname, sep='\t')

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

    # extract the sample name (?)
    temp = df.basename.str.split('_', expand=True)[[0,1]]#.str.join('_')
    df['sample_temp'] = temp[0]+'_'+temp[1]

    # get tech rep numbers -- each mouse has multiple reps
    # and are therefore technical reps
    df['techrep_num'] = df.sort_values(['genotype', 'sample_temp'],
    							ascending=[True, True])\
    							.groupby(['sample_temp']) \
    							.cumcount() + 1

    # get biorep numbers -- each sample is a different mouse
    # and therefore a different biorep
    temp = df[['genotype', 'sample_temp']].drop_duplicates()
    temp.reset_index(inplace=True, drop=True)
    temp['biorep_num'] = temp.sort_values(['genotype', 'sample_temp'],
    							ascending=[True, True])\
    							.groupby(['genotype']) \
    							.cumcount()+1
    df = df.merge(temp, how='left',
                  on=['genotype', 'sample_temp'])

    # sample should be the genotype + mouse id
    # so genotype + biorep
    df['sample'] = df.genotype+'_'+df.biorep_num.astype(str)

    # dataset should be genotype + mouse id + tech rep
    df['dataset'] = df.genotype+'_'+df.biorep_num.astype(str)+'_'+df.techrep_num.astype(str)

    # get the talon run number these will go into
    talon_run_num = 0
    df['talon_run_num'] = np.nan
    for ind, entry, in df.iterrows():
        if ind % datasets_per_run == 0:
            talon_run_num += 1
        df.loc[ind, 'talon_run_num'] = talon_run_num
    df['talon_run_num'] = df.talon_run_num.astype(int)

    return df

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
