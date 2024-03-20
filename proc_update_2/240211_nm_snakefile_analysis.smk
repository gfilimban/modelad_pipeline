import pandas as pd
import os
import sys
import itertools
import swan_vis as swan
import cerberus

from utils import *
from sm_utils import *
from humanized_utils import *

p_dir = '/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update_2/'

configfile: 'config.yml'
config_tsv = 'config.tsv'
p_meta_tsv = 'pseudochromosome_metadata.tsv'
meta_tsv = 'mouse_metadata.tsv'
geno_tsv = 'genotype_metadata.tsv'
an_meta_tsv = '240211_nm_analysis_config.tsv'
auto_dedupe = True

df, p_df = parse_config_file_analysis(config_tsv,
                       meta_tsv,
                       p_meta_tsv,
                       geno_tsv,
                       an_meta_tsv,
                       auto_dedupe=auto_dedupe)

end_modes = ['tss', 'tes']
strands = ['fwd', 'rev']

include: 'cerberus_agg_annot.smk'
include: 'swan.smk'

wildcard_constraints:
    genotype='|'.join([re.escape(x) for x in df.genotype.unique().tolist()]),
    sex='|'.join([re.escape(x) for x in df.sex.unique().tolist()]),
    age='|'.join([re.escape(x) for x in df.age.unique().tolist()]),
    tissue='|'.join([re.escape(x) for x in df.tissue.unique().tolist()]),
    biorep_num='|'.join([re.escape(x) for x in df.biorep_num.unique().tolist()]),
    flowcell='|'.join([re.escape(x) for x in df.flowcell.unique().tolist()]),
    cerberus_run='|'.join([re.escape(x) for x in df.cerberus_run.unique().tolist()]),
    analysis='|'.join([re.escape(x) for x in p_df.analysis.unique().tolist()]),
    feat='|'.join([re.escape(x) for x in ['tss', 'tes', 'ic', 'iso']]),
    obs_col='|'.join([re.escape(x) for x in ['genotype_alias_int', 'genotype_sex']]),
    obs_cond1='|'.join([re.escape(x) for x in p_df.genotype_alias_int.unique().tolist()+p_df.genotype_sex.unique().tolist()]),
    obs_cond2='|'.join([re.escape(x) for x in p_df.genotype_alias_int.unique().tolist()+p_df.genotype_sex.unique().tolist()]),

ruleorder:
    # cerberus_agg_ics_first > cerberus_agg_ics_seq
    # cerb_agg_ics_first > cerb_agg_ics_seq > cerb_agg_ends_first > cerb_agg_ends_seq

rule all:
    input:
        rules.all_swan.input
        # expand(expand(config['analysis']['cerberus']['ca_annot'],
        #         zip,
        #         analysis=p_df.analysis.tolist(),
        #         study=p_df.study.tolist(),
        #         genotype=p_df.genotype.tolist(),
        #         sex=p_df.sex.tolist(),
        #         age=p_df.age.tolist(),
        #         tissue=p_df.tissue.tolist(),
        #         cerberus_run=p_df.cerberus_run.tolist(),
        #         allow_missing=True),
        #         end_mode=['tss', 'tes'])


################################################################################
########################### Cerberus agg + annot ###############################
################################################################################

use rule cerb_agg_ends as cerb_agg_ends_lr with:
    input:
        ref_ends = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['analysis']['cerberus']['agg']['ends'],
                                                  config,
                                                  p_dir),
        ends = p_dir+config['cerberus']['ends']
    params:
        add_ends = True,
        ref = False,
        slack = lambda wc:config['cerberus'][wc.end_mode]['agg_slack'],
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ends = config['analysis']['cerberus']['agg']['ends']

use rule cerberus_agg_ics_cfg as cerb_agg_ics_cfg_lr with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['analysis']['cerberus']['agg']['ics'],
                                                  config,
                                                  p_dir),
        ics = p_dir+config['cerberus']['ics']
    params:
        ref = False,
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        cfg = config['analysis']['cerberus']['agg']['ics_cfg']

use rule cerberus_agg_ics_cli as cerb_agg_ics_cfg_cli_lr with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['analysis']['cerberus']['agg']['ics'],
                                                  config,
                                                  p_dir),
        ics = p_dir+config['cerberus']['ics'],
        cfg = config['analysis']['cerberus']['agg']['ics_cfg']
    output:
        ics = config['analysis']['cerberus']['agg']['ics']

# use rule cerberus_agg_ics as cerb_agg_ics_lr with:
#     input:
#         ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
#                                                   config['analysis']['cerberus']['agg']['ics'],
#                                                   config,
#                                                   p_dir),
#         ics = p_dir+config['cerberus']['ics']
#     params:
#         ref = False,
#         sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
#     output:
#         ics = config['analysis']['cerberus']['agg']['ics']

use rule cerb_write_ref as cerb_write_ref_lr with:
    input:
        ic = config['analysis']['cerberus']['agg']['ics'],
        tss = lambda wc:expand(get_cfg_entries_analysis(wc,
                        p_df,
                        config['analysis']['cerberus']['agg']['ends'])[0],
                        end_mode='tss')[0],
        tes = lambda wc:expand(get_cfg_entries_analysis(wc,
                        p_df,
                        config['analysis']['cerberus']['agg']['ends'])[0],
                        end_mode='tes')[0]
    output:
        h5 = config['analysis']['cerberus']['ca']


use rule cerb_annot as cerb_annot_ref with:
    input:
        h5 = lambda wc:get_final_cerb_entry(wc,
                     p_df,
                     config['analysis']['cerberus']['ca']),
        gtf = p_dir+config['ref']['gtf']
    params:
        source = config['ref']['gtf_ver'],
        gene_source = None
    output:
        h5 = config['analysis']['ref']['cerberus']['ca_annot']

use rule cerb_annot as cerb_annot_run with:
    input:
        h5 = lambda wc:get_final_cerb_entry(wc,
                     p_df,
                     config['analysis']['cerberus']['ca']),
        gtf = p_dir+config['lapa']['filt']['sort_gtf']
    params:
        source = lambda wc:get_df_col(wc, df, 'source'),
        gene_source = None
    output:
        h5 = config['analysis']['cerberus']['ca_annot']

#
use rule cerb_gtf_ids as cerb_update_ref_gtf with:
    input:
        h5 = config['analysis']['ref']['cerberus']['ca_annot'],
        gtf = p_dir+config['ref']['gtf']
    params:
        source = config['ref']['gtf_ver'],
        update_ends = True,
        agg = True
    output:
        gtf = config['analysis']['ref']['cerberus']['gtf']

use rule cerb_gtf_ids as cerb_update_gtf with:
    input:
        h5 = config['analysis']['cerberus']['ca_annot'],
        gtf = p_dir+config['lapa']['filt']['sort_gtf']
    params:
        source = lambda wc:get_df_col(wc, df, 'source'),
        update_ends = True,
        agg = True
    output:
        gtf = config['analysis']['cerberus']['gtf']

use rule cerb_ab_ids as study_cerb_ab with:
    input:
        h5 = config['analysis']['cerberus']['ca_annot'],
        ab = p_dir+config['lapa']['filt']['filt_ab']
    params:
        source = lambda wc:get_df_col(wc, df, 'source'),
        agg = True
    output:
        ab = config['analysis']['cerberus']['ab']
