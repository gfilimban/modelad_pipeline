import pandas as pd
import os
import sys
import swan_vis as swan
import itertools
import cerberus

from utils import *
from sm_utils import *
from humanized_utils import *

configfile: 'config.yml'
config_tsv = 'config.tsv'
p_meta_tsv = 'pseudochromosome_metadata.tsv'
meta_tsv = 'mouse_metadata.tsv'
an_meta_tsv = 'analysis_config.tsv'
auto_dedupe = True

df, p_df = parse_config_file_analysis(config_tsv,
                       meta_tsv,
                       p_meta_tsv,
                       an_meta_tsv,
                       auto_dedupe=auto_dedupe)

end_modes = ['tss', 'tes']
strands = ['fwd', 'rev']

include: 'download.smk'
include: 'samtools.smk'
include: 'refs.smk'
include: 'talon.smk'
include: 'lapa.smk'
include: 'mapping.smk'
include: 'tc.smk'
include: 'cerberus.smk'
include: 'deeptools.smk'
include: 'igvtools.smk'
include: 'pseudochrom.smk'
include: 'swan.smk'

wildcard_constraints:
    genotype='|'.join([re.escape(x) for x in df.genotype.unique().tolist()]),
    sex='|'.join([re.escape(x) for x in df.sex.unique().tolist()]),
    age='|'.join([re.escape(x) for x in df.age.unique().tolist()]),
    tissue='|'.join([re.escape(x) for x in df.tissue.unique().tolist()]),
    biorep_num='|'.join([re.escape(x) for x in df.biorep_num.unique().tolist()]),
    flowcell='|'.join([re.escape(x) for x in df.flowcell.unique().tolist()]),
    cerberus_run='|'.join([re.escape(x) for x in df.cerberus_run.unique().tolist()]),
    obs_col='|'.join([re.escape(x) for x in ['genotype', 'genotype_sex']]),
    obs_cond1='|'.join([re.escape(x) for x in df.genotype.unique().tolist()+df.genotype_sex.unique().tolist()]),
    obs_cond2='|'.join([re.escape(x) for x in df.genotype.unique().tolist()+df.genotype_sex.unique().tolist()]),



ruleorder:
    # cerberus_agg_ics_first > cerberus_agg_ics_seq
    # cerb_agg_ics_first > cerb_agg_ics_seq > cerb_agg_ends_first > cerb_agg_ends_seq

# max_cerb = df.loc[df.cerberus_run.astype(int)==df.cerberus_run.astype(int).max(axis=0)]
# max_cerb = df.loc[df.cerberus_run.astype(int)==3]


rule all:
    input:


################################################################################
################################# Cerberus agg #####################################
################################################################################
################################################################################
################################# Cerberus #####################################
################################################################################

use rule cerb_gtf_to_bed as cerb_get_gtf_ends with:
    input:
        gtf = config['lapa']['filt']['gtf']
    output:
        ends = config['cerberus']['ends']
    params:
        slack = lambda wc:config['cerberus'][wc.end_mode]['slack'],
        dist = lambda wc:config['cerberus'][wc.end_mode]['dist']

use rule cerb_gtf_to_ics as cerb_get_gtf_ics with:
    input:
        gtf = config['lapa']['filt']['gtf']
    output:
        ics = config['cerberus']['ics']

use rule cerb_agg_ends as cerb_agg_ends_lr with:
    input:
        ref_ends = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['cerberus']['agg']['ends'],
                                                  config),
        ends = config['cerberus']['ends']
    params:
        add_ends = True,
        ref = False,
        slack = lambda wc:config['cerberus']['agg'][wc.end_mode]['slack'],
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ends = config['cerberus']['agg']['ends']

use rule cerberus_agg_ics as cerb_agg_ics_lr with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['cerberus']['agg']['ics'],
                                                  config),
        ics = config['cerberus']['ics']
    params:
        ref = False,
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ics = config['cerberus']['agg']['ics']

use rule cerb_write_ref as cerb_write_ref_lr with:
    input:
        ic = config['cerberus']['agg']['ics'],
        tss = lambda wc:expand(get_cfg_entries(wc,
                        p_df,
                        config['cerberus']['agg']['ends']),
                        end_mode='tss')[0],
        tes = lambda wc:expand(get_cfg_entries(wc,
                        p_df,
                        config['cerberus']['agg']['ends']),
                        end_mode='tes')[0]
    output:
        h5 = config['cerberus']['ca']

use rule cerb_annot as cerb_annot_ref with:
    input:
        h5 = config['ref']['cerberus']['ca'],
        gtf = config['ref']['gtf']
    params:
        source = config['ref']['gtf_ver'],
        gene_source = None
    output:
        h5 = config['analysis']['ref']['cerberus']['ca_annot']

use rule cerb_annot as cerb_annot_run with:
    input:
        h5 = config['cerberus']['ca'],
        gtf = config['lapa']['filt']['gtf']
    params:
        source = lambda wc:get_df_col(wc, df, 'source'),
        gene_source = None
    output:
        h5 = config['cerberus']['ca_annot']

#
use rule cerb_gtf_ids as cerb_update_ref_gtf with:
    input:
        h5 = config['analysis']['ref']['cerberus']['ca_annot'],
        gtf = config['ref']['gtf']
    params:
        source = config['ref']['gtf_ver'],
        update_ends = True,
        agg = True
    output:
        gtf = config['analysis']['ref']['cerberus']['gtf']

use rule cerb_gtf_ids as cerb_update_gtf with:
    input:
        h5 = config['cerberus']['ca_annot'],
        gtf = config['lapa']['filt']['gtf']
    params:
        source = lambda wc:get_df_col(wc, df, 'source'),
        update_ends = True,
        agg = True
    output:
        gtf = config['cerberus']['gtf']

use rule cerb_ab_ids as study_cerb_ab with:
    input:
        h5 = config['cerberus']['ca_annot'],
        ab = config['lapa']['filt']['filt_ab']
    params:
        source = lambda wc:get_df_col(wc, df, 'source'),
        agg = True
    output:
        ab = config['cerberus']['ab']
