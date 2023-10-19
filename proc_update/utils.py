import pysam
import pandas as pd
import re
import math
import gzip
import pyranges as pr
import numpy as np
import scanpy as sc
import swan_vis as swan

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


def ab_add_rescue_ism_cat(ab):
    """
    Update LAPA abundance w/ ISM rescue category for those
    that were assigned new ends from LAPA
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
    df = pr.read_gtf(gtf, as_df=True)
    df = df.loc[~df.Chromosome.str.contains('SIRV')]
    df = df.loc[~df.Chromosome.str.contains('ERCC')]
    df = df.loc[df.Feature == 'transcript']
    filt_df = df[['gene_id', 'transcript_id']].drop_duplicates()
    filt_df = filt_df.rename({'gene_id':'gid',
                              'transcript_id':'tid'}, axis=1)
    return filt_df

def filter_lapa(ab,
                gtf,
                t_novs,
                g_novs,
                filt_spikes,
                ofile):
    """
    Filter LAPA transcripts and output a pass list of
    passed transcripts that are in the list of valid
    transcript novelites, gene novelties, and depending
    on whether or not we're filtering out spike-ins.

    Parameters:
        ab (str): Path to abundance file
        gtf (str): Path to GTF
        t_novs (list of str): List of transcript novelties
            to retain transcripts from
        g_novs (list of str): List of gene novelties to
            retain transcripts from
        filt_spikes (bool): Whether to remove spikeins
        ofile (str): Path / name of output file
    """

    # filter based on novelty after defining
    # rescue ISMS
    df = ab_add_rescue_ism_cat(ab)
    filt_df = filter_lapa_on_nov(df,
                                 t_novs,
                                 g_novs)

    # filter out spike-ins
    if filt_spikes:
        temp = filter_spikes(gtf)
        filt_df = filt_df.merge(temp, how='inner')

    filt_df.to_csv(ofile, index=False, sep='\t')
