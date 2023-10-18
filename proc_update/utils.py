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
    """
    Process metadat information. Remove spaces and forward
    slashes with underscores. Verify that mouse ids aren't
    duplicated.
    """
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

    ############ Basename + fname
    df['basename'] = df.fname.str.rsplit('/', n=1, expand=True)[1]
    df['path'] = df.fname.str.rsplit('/', n=1, expand=True)[0]

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

    df['flowcell'] = df.flowcell.astype(str)
    df['biorep_num'] = df.biorep_num.astype(str)

    return df

def reverse_alignment(infile, outfile, threads=1):
    """
    Flip alignents that are in the reverse orientation
    (ie are oriented 3'->5'), such that they are in the
    forward orientation (5'->3')
    """

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
