import pandas as pd
import os
import sys
import swan_vis as swan
import itertools
import cerberus

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from utils import *

# settings we can change each time it's run
configfile: 'config.yml'
config_tsv = '230607_config.tsv'
cerb_tsv = 'cerberus.tsv' # gtf_to_bed and agg_ends settings
datasets_per_run = 4 # number of datasets per talon run
auto_dedupe = True # deduplicate runs w/ same stem but different chop numbers
cerb_settings = pd.read_csv(cerb_tsv, sep='\t')

# should be static
meta_tsv = 'mouse_metadata.tsv'


df, dataset_df = parse_config_file(config_tsv,
                       meta_tsv,
                       datasets_per_run=datasets_per_run,
                       auto_dedupe=auto_dedupe)

# from the df w/ one line for each file
batches = df.batch.tolist()
datasets = df.dataset.tolist()
flowcells = df.flowcell.tolist()
samples = df['sample'].tolist()
platforms = df.platform.tolist()
genotypes = df.genotype.unique().tolist()

# from the df w/ one line for each dataset / mouse
temp = dataset_df[['study', 'talon_run_num']].drop_duplicates()
temp = temp.groupby('study').max().reset_index()
max_talon_runs = temp.talon_run_num.tolist()
studies = temp.study.tolist()

# for cerberus annotate_transcriptome
sources = ['vM21']+studies
source_df = pd.DataFrame()
source_df['source'] = sources

end_modes = ['tss', 'tes']

sr_files = ['/share/crsp/lab/model-ad/share/bulkRNA/5x_GWAS/5xBin1/5xBin1_4mo/hipp_F_BIN1_HO_4mo_13047_S16-STARAligned.out.bam',
            '/share/crsp/lab/model-ad/share/bulkRNA/5x_GWAS/5xBin1/5xBin1_4mo/hipp_F_5xFADHEMI_BIN1_HO_4mo_13019_S1-STARAligned.out.bam',
            '/share/crsp/lab/model-ad/share/bulkRNA/5x_GWAS/5xGWAScontrols_4mo_2021July/hipp_F_5xFADWT_4mo_11627_S43-STARAligned.out.bam',
            '/share/crsp/lab/model-ad/share/bulkRNA/5x_GWAS/5xGWAScontrols_4mo_2021July/hipp_F_5xFADHEMI_4mo_11616_S21-STARAligned.out.bam']
sr_df = pd.DataFrame(data=sr_files, columns=['fname'])
sr_df['mouse_id'] = sr_df['fname'].str.rsplit('_', n=3, expand=True)[2]
mouse_ids = sr_df.mouse_id.tolist()

# lr
lr_datasets = ['5xBIN1_HO_F_4_months_HC_1',
               'BIN1_HO_F_4_months_HC_2',
               '5xFADHEMI_F_4_months_HC_1',
               '5xFADWT_F_4_months_HC_1']

strands = ['fwd', 'rev']

def get_genotype_pairs(df, pair_num):
    genotypes = df.genotype.unique().tolist()
    pairs = list(itertools.combinations(genotypes, 2))
    if pair_num == 0:
        g = [p[0] for p in pairs]
    elif pair_num == 1:
        g = [p[1] for p in pairs]
    return g

def get_df_col(wc, df, col):
    val = df.loc[(df.dataset==wc.dataset)&(df.flowcell==wc.flowcell), col].values[0]
    return val

def get_dataset_df_col(wc, df, col):
    val = df.loc[df.dataset==wc.dataset, col].values[0]
    return val

# config formatting errors
if len(df.batch.unique()) > 1:
    raise ValueError('Must only have one batch per config')
else:
    batch = batches[0]

wildcard_constraints:
    genotype1= '|'.join([re.escape(x) for x in genotypes]),
    genotype2= '|'.join([re.escape(x) for x in genotypes]),
    batch=batch,
    mouse_ids= '|'.join([re.escape(x) for x in mouse_ids])

ruleorder:
    first_talon > seq_talon

rule all:
    input:
        # curr working
        # expand(config['data']['map_stats'],
        #    zip,
        #    batch=batches,
        #    dataset=datasets,
        #    flowcell=flowcells),
        # expand(config['data']['tc_stats'],
        #   zip,
        #   batch=batches,
        #   dataset=datasets,
        #   flowcell=flowcells),
        # expand(expand(config['data']['lapa_filt_ab'],
        #        zip,
        #        study=studies,
        #        allow_missing=True),
        #        batch=batch),
        # expand(expand(config['data']['ca_annot_2'],
        #        zip,
        #        study=studies,
        #        allow_missing=True),
        #        batch=batch),
        # expand(expand(config['data']['sg'],
        #        zip,
        #        study=studies,
        #        allow_missing=True),
        #        batch=batch),
        # expand(config['sr']['bw'],
        #        mouse_id=mouse_ids,
        #        strand=strands),
        # expand(config['data']['bw'],
        #        batch=batch,
        #        dataset=lr_datasets,
        #        strand=strands)
        # expand(config['data']['ca_ref_gtf'],
        #        zip,
        #        batch=batch),
        # expand(expand(config['data']['ca_gtf'],
        #        zip,
        #        study=studies,
        #        allow_missing=True),
        #        batch=batch),
        # expand(expand(config['data']['ca_ab'],
        #        zip,
        #        study=studies,
        #        allow_missing=True),
        #        batch=batch),

        # expand(expand(config['data']['ca_annot_2'],
        #        zip,
        #        study=studies,
        #        allow_missing=True),
        #        batch=batch),

        # trying to figure out stupid sequential Cerberus
        # expand(config['data']['ca_annot'],
        #        zip,
        #        batch=batch,
        #        study=source_df.loc[source_df.index.max(), 'source'],
        #        cerb_run=source_df.index.max()),
        # expand(config['data']['ca_annot'],
        #        zip,
        #        batch=batch,
        #        study=source_df.loc[1, 'source'],
        #        cerb_run=1),
        # expand(config['data']['ca_annot'],
        #        zip,
        #        batch=batch,
        #        study=source_df.loc[0, 'source'],
        #        cerb_run=0),

        # need to clean up these guyes
        # expand(config['data']['sg'], batch=batches),
        # expand(expand(config['data']['die_tsv'],
        #        zip,
        #        genotype1=get_genotype_pairs(df, 0),
        #        genotype2=get_genotype_pairs(df, 1),
        #        allow_missing=True),
        #        batch=batches[0]),
        # expand(expand(config['data']['de_tsv'],
        #        zip,
        #        genotype1=get_genotype_pairs(df, 0),
        #        genotype2=get_genotype_pairs(df, 1),
        #        allow_missing=True),
        #        batch=batches[0],
        #        feature=['gene', 'iso'])
        # expand(expand(config['data']['de_tsv'],
        #        zip,
        #        genotype1=get_genotype_pairs(df, 0)[0],
        #        genotype2=get_genotype_pairs(df, 1)[0],
        #        allow_missing=True),
        #        batch=batches[0],
        #        feature=['gene', 'iso'])
        expand(expand(config['data']['bam_label_merge'],
               zip,
               dataset=datasets,
               allow_missing=True),
               batch=batch)

################################################################################
########################### Ref. processing ####################################
################################################################################

rule dl:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget -O {output.out} {params.link}"

rule gunzip:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "gunzip -c {input.gz} > {output.out}"

use rule dl as dl_gencode_gtf with:
  params:
    link = config['ref']['gtf_link']
  output:
    out = config['ref']['gtf_gz']

use rule gunzip as gunzip_gencode with:
  input:
    gz = config['ref']['gtf_gz']
  output:
    out = config['ref']['gtf']

use rule dl as dl_ref with:
    params:
        link = config['ref']['fa_link']
    output:
        out = config['ref']['fa_gz']

use rule gunzip as gunzip_ref with:
    input:
        gz = config['ref']['fa_gz']
    output:
        out = config['ref']['fa']

rule get_chrom_sizes:
    input:
        fa = config['ref']['fa']
    resources:
        threads = 1,
        mem_gb = 8
    output:
        chrom_sizes = config['ref']['chrom_sizes']
    shell:
        """
        faidx {input.fa} -i chromsizes > {output.chrom_sizes}
        """

rule get_utr_fix_gtf:
    input:
        gtf = config['ref']['gtf']
    resources:
        threads = 1,
        mem_gb = 8
    output:
        gtf = config['ref']['gtf_utr']
    shell:
        """
        gencode_utr_fix \
            --input_gtf {input.gtf} \
            --output_gtf {output.gtf}
        """
rule get_annot_sjs:
    input:
        gtf = config['ref']['gtf'],
        fa = config['ref']['fa']
    resources:
        threads = 4,
        mem_gb = 16
    params:
        tc = config['tc_path']
    output:
        sjs = config['ref']['sjs']
    shell:
        """
        python {params.tc}accessory_scripts/get_SJs_from_gtf.py \
             --f {input.gtf} \
             --g {input.fa} \
             --minIntronSize 21 \
             --o {output.sjs}
        """

use rule dl as dl_ca with:
    params:
        link = config['ref']['ca_link']
    output:
        out = config['ref']['ca']

rule ca_to_tsv:
    input:
        ca = config['ref']['ca']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        opref = config['ref']['ca_ics'].split('_ic.tsv')[0]
    output:
        expand(config['ref']['ca_ends'],
               zip,
               end_mode=end_modes),
        config['ref']['ca_ics']
    run:
        cerberus.write_h5_to_tsv(input.ca,
                                 params.opref)

################################################################################
############################# General rules ####################################
################################################################################
rule alignment_stats:
    resources:
        threads = 2,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools stats {input.alignment} | grep ^SN | cut -f 2- | grep -e 'reads map
ped' -e 'reads unmapped' -e 'average length' -e 'maximum length' | sed '/reads mapped and paired/d' > {output.stats}
        """

################################################################################
################################ Mapping #######################################
################################################################################
# rule rev_fastq:
#     input:
#         fastq = lambda wc: get_df_col(wc, df, 'fname')
#     resources:
#         threads = 1,
#         mem_gb = 16
#     params:
#     output:
#         fastq = config['data']['fastq']
#     run:
#         flip_fastq(input.fastq, output.fastq)

rule map:
  resources:
    threads = 32,
    mem_gb = 64
  shell:
      """
      module load minimap2
      minimap2 --MD \
     			-t {resources.threads} \
     			-ax splice \
     			-k14 \
                --junc-bed {input.sjs} \
     		    {input.ref_fa} {input.fastq} > {output.sam} 2> {output.log}
    """

rule sam_to_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools sort \
            --threads {resources.threads} \
            -O bam {input.sam} > {output.bam}
        samtools index -@ {resources.threads} {output.bam}
        """

rule merge_alignment:
    resources:
        threads = 32,
        mem_gb = 64
    shell:
        """
        module load samtools
        samtools merge -o {output.bam} {input.files}
        """

use rule map as map_reads with:
    input:
        fastq = lambda wc: get_df_col(wc, df, 'fname'),
        ref_fa = config['ref']['fa'],
        sjs = config['ref']['sjs']
    output:
        sam = temporary(config['data']['sam']),
        log = config['data']['sam_log']

# use rule sam_to_bam as convert_sam:
#     input:
#         sam = config['data']['sam']
#     output:
#         bam = config['data']['bam']

use rule alignment_stats as map_stats with:
    input:
        alignment = config['data']['sam']
    output:
        stats = config['data']['map_stats']

rule rev_alignment:
    input:
        sam = config['data']['sam']
    resources:
        threads = 8,
        mem_gb = 32
    output:
        sam_rev = temporary(config['data']['sam_rev'])
    run:
        reverse_alignment(input.sam, output.sam_rev, resources.threads)

################################################################################
############################ TranscriptClean ###################################
################################################################################
rule tc:
    resources:
        mem_gb = 80,
        threads = 16
    shell:
        """
        python {params.tc}TranscriptClean.py \
            -t {resources.threads} \
            --sam {input.sam} \
            --genome {input.fa} \
            --canonOnly \
            --primaryOnly \
            --deleteTmp \
            --correctMismatches False \
            --correctIndels True \
            --tmpDir {params.opref}_temp/ \
            --outprefix {params.opref}
        """

use rule tc as tc_sam with:
    input:
        sam = config['data']['sam_rev'],
        fa = config['ref']['fa']
    params:
        tc = config['tc_path'],
        opref = config['data']['sam_clean'].rsplit('_clean.sam', maxsplit=1)[0]
    output:
        sam = temporary(config['data']['sam_clean'])

use rule alignment_stats as tc_stats with:
    input:
        alignment = config['data']['sam_clean']
    output:
        stats = config['data']['tc_stats']

################################################################################
############################## TALON label #####################################
################################################################################
rule talon_label:
    input:
        fa = config['ref']['fa'],
        sam = config['data']['sam_clean']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        opref = config['data']['sam_label'].rsplit('_labeled.sam', maxsplit=1)[0]
    output:
        sam = temporary(config['data']['sam_label'])
    shell:
        """
        talon_label_reads \
            --f {input.sam} \
            --g {input.fa} \
            --tmpDir {params.opref} \
            --ar 20 \
            --deleteTmp \
            --o {params.opref}
        """

use rule sam_to_bam as talon_label_bam with:
    input:
        sam = config['data']['sam_label']
    output:
        bam = config['data']['bam_label_sorted']

################################################################################
################################# TALON ########################################
################################################################################

# merge files from the different flowcell into the same talon input file
def get_merge_talon_label_files(wc, batches, df, config_entry):
    temp = df.loc[df.dataset == wc.dataset]
    datasets = temp.dataset.tolist()
    flowcells = temp.flowcell.tolist()
    files = expand(config_entry,
                   zip,
                   batch=batches,
                   dataset=datasets,
                   flowcell=flowcells)
    return files

use rule merge_alignment as merge_talon_label with:
    input:
        files = lambda wc:get_merge_talon_label_files(wc,
                                                      batches,
                                                      df,
                                                      config['data']['bam_label_sorted'])
    output:
        bam = config['data']['bam_label_merge']

def get_talon_run_info(wc, df, config_entry=None, col=False):
    temp = df.loc[df.talon_run_num == int(wc.talon_run)]
    datasets = get_study_datasets(wc, temp)
    temp = temp.loc[temp.dataset.isin(datasets)]
    if config_entry:
        files = expand(expand(config_entry,
                       zip,
                       dataset=datasets,
                       allow_missing=True),
                       batch=batch)
        temp['file'] = files
    if col:
        return temp[col].tolist()
    else:
        return temp

rule talon_config:
    input:
        files = lambda wc:get_talon_run_info(wc, dataset_df, config['data']['bam_label_merge'], col='file')
    resources:
        threads = 1,
        mem_gb = 1
    params:
        df = dataset_df
    output:
        config = config['data']['talon_config']
    run:
        config = get_talon_run_info(wildcards, params.df)
        config = config[['dataset', 'sample', 'platform']].copy(deep=True)
        config['file'] = input.files
        config.to_csv(output.config, header=None, sep=',', index=False)

rule talon_init:
	input:
		ref_gtf = config['ref']['gtf']
	output:
		db = config['ref']['talon_db']
	params:
		talon_opref = config['ref']['talon_db'].rsplit('.db', maxsplit=1)[0],
		genome = 'mm10',
		annot = 'vM21'
	resources:
		mem_gb = 32,
		threads = 16
	shell:
		'talon_initialize_database \
    		--f {input.ref_gtf} \
    		--g {params.genome} \
    		--a {params.annot} \
    		--l 0 \
    		--idprefix TALON \
    		--5p 500 \
    		--3p 300 \
    		--o {params.talon_opref}'

rule talon:
    resources:
        mem_gb = 256,
        threads = 30
    shell:
        """
        ref_db={input.ref}_{wildcards.batch}_{wildcards.study}
        cp {input.ref} ${{ref_db}}
        talon \
            --f {input.config} \
            --db ${{ref_db}} \
            --build {params.genome} \
            --tmpDir {params.opref}_temp/ \
            --threads {resources.threads} \
            --o {params.opref}
        mv ${{ref_db}} {params.opref}_talon.db
        """

use rule talon as first_talon with:
    input:
        ref = config['ref']['talon_db'],
        config = expand(config['data']['talon_config'],
                        zip,
                        talon_run=1,
                        allow_missing=True)[0]
    params:
        genome = 'mm10',
        opref = expand(config['data']['talon_db'],
                      zip,
                      talon_run=1,
                      allow_missing=True)[0].rsplit('_talon', maxsplit=1)[0],
    output:
        db = expand(config['data']['talon_db'],
                   zip,
                   talon_run=1,
                   allow_missing=True)[0],
        annot = expand(config['data']['read_annot'],
                                  zip,
                                  talon_run=1,
                                  allow_missing=True)[0]

use rule talon as seq_talon with:
    input:
        ref = lambda wc: expand(config['data']['talon_db'],
                        batch=wc.batch,
                        study=wc.study,
                        talon_run=int(wc.talon_run)-1)[0],
        config = config['data']['talon_config']
    params:
        genome = 'mm10',
        opref = config['data']['talon_db'].rsplit('_talon', maxsplit=1)[0]
    output:
        db = config['data']['talon_db'],
        annot = config['data']['read_annot']

def get_talon_max_output_file(wc, df, config_entry):
    datasets = get_study_datasets(wc, df)
    temp = df.loc[df.dataset.isin(datasets)]
    n = temp.talon_run_num.max()
    file = expand(config_entry,
                  zip,
                  batch=wc.batch,
                  study=wc.study,
                  talon_run=n)[0]
    return file

rule talon_read_annot:
    input:
        db = lambda wc: get_talon_max_output_file(wc, dataset_df, config['data']['talon_db'])
    resources:
        mem_gb = 32,
        threads = 1
    params:
        genome = 'mm10',
        opref = config['data']['full_read_annot'].rsplit('_talon', maxsplit=1)[0]
    output:
        tsv = config['data']['full_read_annot']
    shell:
        """
        talon_fetch_reads \
            --db {input.db} \
            --build {params.genome} \
            --o {params.opref}
        """

rule talon_unfilt_ab:
    input:
        db = lambda wc: get_talon_max_output_file(wc,
                                                  dataset_df,
                                                  config['data']['talon_db'])
    resources:
        threads = 1,
        mem_gb = 32
    params:
        genome = 'mm10',
        annot = 'vM21',
        opref = config['data']['ab'].rsplit('_talon', maxsplit=1)[0]
    output:
        ab = config['data']['ab']
    shell:
        """
        talon_abundance \
            --db {input.db} \
            -a {params.annot} \
            -b {params.genome} \
            --o {params.opref}
        """

rule talon_filt:
    input:
        db = lambda wc: get_talon_max_output_file(wc,
                                                  dataset_df,
                                                  config['data']['talon_db'])
    resources:
        threads = 1,
        mem_gb = 128
    params:
        annot = 'vM21'
    output:
        list = config['data']['filt_list']
    shell:
        """
        talon_filter_transcripts \
            --db {input.db} \
            -a {params.annot} \
            --maxFracA=0.5 \
            --minCount=5 \
            --minDatasets=2 \
            --o {output.list}
        """

rule talon_filt_ab:
    input:
        db = lambda wc: get_talon_max_output_file(wc,
                                                  dataset_df,
                                                  config['data']['talon_db']),
        filt = config['data']['filt_list']
    resources:
        threads = 1,
        mem_gb = 128
    params:
        genome = 'mm10',
        annot = 'vM21',
        opref = config['data']['filt_ab'].rsplit('_talon', maxsplit=1)[0]
    output:
        filt_ab = config['data']['filt_ab']
    shell:
        """
        talon_abundance \
            --db {input.db} \
            -a {params.annot} \
            -b {params.genome} \
            --whitelist {input.filt} \
            --o {params.opref}
        """

rule talon_gtf:
    input:
        db = lambda wc: get_talon_max_output_file(wc,
                                                  dataset_df,
                                                  config['data']['talon_db']),
        filt = config['data']['filt_list']
    resources:
        threads = 1,
        mem_gb = 128
    params:
        genome = 'mm10',
        annot = 'vM21',
        opref = config['data']['filt_gtf'].rsplit('_talon', maxsplit=1)[0]
    output:
        gtf = config['data']['filt_gtf']
    shell:
        """
        talon_create_GTF \
            --db {input.db} \
            -b {params.genome} \
            -a {params.annot} \
            --whitelist {input.filt} \
            --observed \
            --o {params.opref}
        """

################################################################################
################################# LAPA #########################################
################################################################################
def get_lapa_settings(wc, lapa_ends, kind):
    lapa_ends = expand(lapa_ends,
                       zip,
                       batch=wc.batch,
                       study=wc.study,
                       end_mode=wc.end_mode)[0]
    if kind == 'temp_file':
        temp = os.path.dirname(lapa_ends)+'/'
        if wc.end_mode == 'tes':
            temp += 'polyA_clusters.bed'
        elif wc.end_mode == 'tss':
            temp += 'tss_clusters.bed'
        return temp
    elif kind == 'lapa_cmd':
        if wc.end_mode == 'tes':
            return 'lapa'
        elif wc.end_mode == 'tss':
            return 'lapa_tss'

def get_study_datasets(wc, dataset_df):
    datasets = dataset_df.loc[dataset_df.study==wc.study, 'dataset'].tolist()
    return datasets

def get_lapa_run_info(wc, df, config_entry=None, col=False):
    temp = df.copy(deep=True)
    datasets = get_study_datasets(wc, temp)
    temp = temp.loc[temp.dataset.isin(datasets)]
    if config_entry:
        files = expand(expand(config_entry,
                       zip,
                       dataset=datasets,
                       allow_missing=True),
                       batch=batch)
        temp['file'] = files
    if col:
        return temp[col].tolist()
    else:
        return temp

rule lapa_config:
    input:
        files = lambda wc: get_lapa_run_info(wc, dataset_df, config['data']['bam_label_merge'], 'file')
    resources:
        threads = 1,
        mem_gb = 1
    params:
        df = dataset_df
    output:
        config = config['data']['lapa_config']
    run:
        config = get_lapa_run_info(wildcards, params.df)
        config = config[['sample', 'dataset']].copy(deep=True)
        config['file'] = input.files
        config.columns = ['sample', 'dataset', 'path']
        config.to_csv(output.config, sep=',', index=False)

rule lapa_call_ends:
    input:
        config = config['data']['lapa_config'],
        fa = config['ref']['fa'],
        gtf = config['ref']['gtf_utr'],
        chrom_sizes = config['ref']['chrom_sizes']
    resources:
        threads = 4,
        mem_gb = 32
    params:
        opref = config['data']['lapa_ends'].rsplit('/', maxsplit=1)[0]+'/',
        lapa_cmd = lambda wc: get_lapa_settings(wc, config['data']['lapa_ends'], 'lapa_cmd'),
        lapa_end_temp = lambda wc: get_lapa_settings(wc, config['data']['lapa_ends'], 'temp_file'),
    output:
        ends = config['data']['lapa_ends']
    shell:
        """
        rm -rf {params.opref}
        {params.lapa_cmd} \
            --alignment {input.config} \
            --fasta {input.fa} \
            --annotation {input.gtf} \
            --chrom_sizes {input.chrom_sizes} \
            --output_dir {params.opref}
        if [ {params.lapa_end_temp} != {output.ends} ]
        then
            cp {params.lapa_end_temp} {output.ends}
        fi
        """

rule lapa_link:
    input:
        annot = config['data']['full_read_annot'],
        tss = expand(config['data']['lapa_ends'], end_mode='tss', allow_missing=True)[0],
        tes = expand(config['data']['lapa_ends'], end_mode='tes', allow_missing=True)[0]
    resources:
        threads = 1,
        mem_gb = 256
    params:
        tss_dir = expand(config['data']['lapa_ends'], end_mode='tss', allow_missing=True)[0].rsplit('/', maxsplit=1)[0]+'/',
        tes_dir = expand(config['data']['lapa_ends'], end_mode='tes', allow_missing=True)[0].rsplit('/', maxsplit=1)[0]+'/'
    output:
        links = config['data']['lapa_links']
    shell:
        """
        lapa_link_tss_to_tes \
            --alignment {input.annot} \
            --lapa_dir {params.tes_dir} \
            --lapa_tss_dir {params.tss_dir} \
            --output {output.links}
        """

rule lapa_correct_talon:
    input:
        gtf = config['data']['filt_gtf'],
        ab = config['data']['filt_ab'],
        annot = config['data']['full_read_annot'],
        links = config['data']['lapa_links']
    resources:
        threads = 1,
        mem_gb = 64
    output:
        gtf = config['data']['lapa_gtf'],
        ab = config['data']['lapa_ab']
    shell:
        """
        lapa_correct_talon \
                --links {input.links} \
                --read_annot {input.annot} \
                --gtf_input {input.gtf} \
                --gtf_output {output.gtf} \
                --abundance_input {input.ab} \
                --abundance_output {output.ab} \
                --keep_unsupported
        """

rule filt_lapa:
    input:
        ab = config['data']['lapa_ab'],
        gtf = config['data']['lapa_gtf']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        t_novs = ['Known', 'NIC', 'NNC', 'ISM_rescue'],
        g_novs = ['Known'],
        filt_spikes = True
    output:
        filt_list = config['data']['lapa_filt_list']
    run:
        # filter based on novelty after defining
        # rescue ISMS
        df = ab_add_rescue_ism_cat(input.ab)
        filt_df = filter_lapa_on_nov(df,
                                     params.t_novs,
                                     params.g_novs)

        # filter out spike-ins
        if params.filt_spikes:
            temp = filter_spikes(input.gtf)
            filt_df = filt_df.merge(temp, how='inner')

        filt_df.to_csv(output.filt_list, index=False, sep='\t')

rule filt_lapa_ab:
    input:
        ab = config['data']['lapa_ab'],
        filt_list = config['data']['lapa_filt_list']
    resources:
        threads = 4,
        mem_gb = 32
    output:
        ab = config['data']['lapa_filt_ab']
    run:
        df = filt_lapa_ab(input.ab,
                          input.filt_list)
        df.to_csv(output.ab, sep='\t', index=False)



rule filt_lapa_gtf:
    input:
        gtf = config['data']['lapa_gtf'],
        filt_list = config['data']['lapa_filt_list']
    resources:
        threads = 4,
        mem_gb = 32
    output:
        gtf = config['data']['lapa_filt_gtf']
    run:
        gtf = filt_lapa_gtf(input.gtf,
                            input.filt_list)
        gtf.to_gtf(output.gtf)


################################################################################
############################## Cerberus ########################################
################################################################################
def get_cerb_settings(wc, df, col):
    """
    Get Cerberus end-related settings
    """
    col = f'{wc.end_mode}_{col}'
    return df[col].tolist()[0]

rule cerb_gtf_to_bed:
    input:
        gtf = config['data']['lapa_filt_gtf']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        ends = config['data']['ends']
    params:
        slack = lambda wc:get_cerb_settings(wc, cerb_settings, 'slack'),
        dist = lambda wc:get_cerb_settings(wc, cerb_settings, 'dist')
    run:
        cerberus.gtf_to_bed(input.gtf,
                            wildcards.end_mode,
                            output.ends,
                            dist=params.dist,
                            slack=params.slack)

rule cerb_gtf_to_ics:
    input:
        gtf = config['data']['lapa_filt_gtf']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        ics = config['data']['ics']
    run:
        cerberus.gtf_to_ics(input.gtf,
                            output.ics)


def get_agg_settings(wc, param='files'):
    settings = {}
    files = []
    sources = []

    if 'end_mode' in wc.keys():
        ic = False
    else:
        ic = True

    # get reference ics / tsss / tess first
    # first source is 'cerberus' to indicate
    # preservation of source / novelty
    if ic == True:
        files += [config['ref']['ca_ics']]
    else:
        files += expand(config['ref']['ca_ends'],
                        zip,
                        end_mode=wc.end_mode)
    sources += ['cerberus']

    # then get the files from each study
    # use the studies as the sources
    if ic == True:
        files += expand(expand(config['data']['ics'],
               zip,
               study=studies,
               allow_missing=True),
               batch=batch)
    else:
        files += expand(expand(config['data']['ends'],
               zip,
               study=studies,
               allow_missing=True),
               batch=batch,
               end_mode=wc.end_mode)
    sources += studies

    settings['file'] = files
    settings['source'] = sources

    return settings[param]

rule cerb_agg_ends:
    input:
        files = lambda wc:get_agg_settings(wc, 'file')
    resources:
      threads = 4,
      mem_gb = 32
    params:
        add_ends = True,
        refs = False,
        slack = lambda wc:get_cerb_settings(wc, cerb_settings, 'agg_slack'),
        sources = lambda wc:get_agg_settings(wc, 'source')
    output:
        ends = config['data']['agg_ends']
    run:
        refs = [params.refs for i in range(len(input.files))]
        add_ends = [params.add_ends for i in range(len(input.files))]
        cerberus.agg_ends(input.files,
                          add_ends,
                          refs,
                          params.sources,
                          wildcards.end_mode,
                          params.slack,
                          output.ends)

rule cerb_agg_ics:
  input:
      files = lambda wc:get_agg_settings(wc, 'file')
  resources:
    threads = 4,
    mem_gb = 32
  params:
      refs = False,
      sources = lambda wc:get_agg_settings(wc, 'source')
  output:
      ics = config['data']['agg_ics']
  run:
      refs = [params.refs for i in range(len(input.files))]
      cerberus.agg_ics(input.files,
                        refs,
                        params.sources,
                        output.ics)

rule cerb_write_ref:
    input:
        ic = config['data']['agg_ics'],
        tss = lambda wc:expand(config['data']['agg_ends'],
                               batch=wc.batch,
                               end_mode='tss')[0],
        tes = lambda wc:expand(config['data']['agg_ends'],
                              batch=wc.batch,
                              end_mode='tes')[0]
    resources:
        threads = 4,
        mem_gb = 64
    output:
        h5 = config['data']['ca_ref']
    run:
        cerberus.write_reference(input.tss,
                                 input.tes,
                                 input.ic,
                                 output.h5)


################################################################################
######################### Cerberus annot + ID replacement ######################
################################################################################

rule cerb_annot:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.annotate_transcriptome(input.gtf,
                                        input.h5,
                                        params.source,
                                        params.gene_source,
                                        output.h5)

use rule cerb_annot as ref_cerb_annot with:
    input:
        h5 = config['data']['ca_ref'],
        gtf = config['ref']['gtf']
    params:
        source = 'vM21',
        gene_source = None
    output:
        h5 = config['data']['ca_ref_annot']

use rule cerb_annot as study_cerb_annot with:
    input:
        h5 = config['data']['ca_ref_annot'],
        gtf = config['data']['lapa_filt_gtf']
    params:
        source = lambda wc:wc.study,
        gene_source = 'vM21'
    output:
        h5 = config['data']['ca_annot_2']

# # first one should be with vM21
# use rule cerb_annot as first_cerb_annot with:
#     input:
#         h5 = expand(config['data']['ca_ref'],
#                     batch=batch)[0],
#         gtf = config['ref']['gtf']
#     params:
#         source = 'vM21',
#         gene_source = None
#     output:
#         h5 = expand(config['data']['ca_annot'],
#                zip,
#                batch=batch,
#                study=source_df.loc[0, 'source'],
#                cerb_run=0)[0]
#
# def get_seq_cerb_h5(wc, df):
#     cerb_run = int(wc.cerb_run)-1
#     if cerb_run < 0:
#         cerb_run = 0
#     study = df.loc[cerb_run, 'source']
#     h5 = expand(config['data']['ca_annot'],
#                            zip,
#                            batch=wc.batch,
#                            study=study,
#                            cerb_run=cerb_run)[0]
#     return h5
#
# use rule cerb_annot as seq_cerb_annot with:
#     input:
#         h5 = lambda wc: get_seq_cerb_h5(wc, source_df),
#         gtf = config['data']['lapa_filt_gtf']
#     params:
#         source = lambda wc:wc.study,
#         gene_source = 'vM21'
#     output:
#         h5 = config['data']['ca_annot']


rule cerb_gtf_ids:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.replace_gtf_ids(input.h5,
                                 input.gtf,
                                 params.source,
                                 True,
                                 True,
                                 output.gtf)

use rule cerb_gtf_ids as ref_cerb_gtf with:
    input:
        h5 = config['data']['ca_ref_annot'],
        gtf = config['ref']['gtf']
    params:
        source = 'vM21',
    output:
        gtf = config['data']['ca_ref_gtf']

use rule cerb_gtf_ids as study_cerb_gtf with:
    input:
        h5 = config['data']['ca_annot_2'],
        gtf = config['data']['lapa_filt_gtf']
    params:
        source = lambda wc:wc.study,
    output:
        gtf = config['data']['ca_gtf']

rule cerb_ab_ids:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.replace_ab_ids(input.ab,
                                input.h5,
                                params.source,
                                True,
                                output.ab)

use rule cerb_ab_ids as study_cerb_ab with:
    input:
        h5 = config['data']['ca_annot_2'],
        ab = config['data']['lapa_filt_ab']
    params:
        source = lambda wc:wc.study,
    output:
        ab = config['data']['ca_ab']

################################################################################
################################ Swan ##########################################
################################################################################

rule save_swan_metadata:
    resources:
        mem_gb = 8,
        threads = 1
    output:
        meta = config['data']['swan_meta']
    run:
        dataset_df.to_csv(output.meta, sep='\t', index=False)


def make_sg(input, params, wildcards):

    # initialize
    sg = swan.SwanGraph()
    sg.add_annotation(input.annot)
    sg.add_transcriptome(input.gtf, include_isms=True)
    sg.add_abundance(input.filt_ab)
    sg.add_abundance(input.ab, how='gene')
    sg.save_graph(params.prefix)

    sg.add_metadata(input.meta)
    # add metadata
    # adatas = sg.get_adatas()
    # for adata in adatas:
    #     adata.obs['genotype'] = adata.obs.dataset.str.rsplit('_', n=2, expand=True)[0]

    sg.save_graph(params.prefix)

rule make_sg:
    input:
        annot = config['data']['ca_ref_gtf'],
        gtf = config['data']['ca_gtf'],
        filt_ab = config['data']['ca_ab'],
        ab = config['data']['ab'],
        meta = config['data']['swan_meta']
    params:
        prefix = config['data']['sg'].replace('.p', '')
    resources:
        mem_gb = 128,
        threads = 1
    output:
        sg = config['data']['sg']
    run:
        make_sg(input, params, wildcards)

rule swan_die:
    input:
        sg = config['data']['sg']
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = config['data']['die_tsv']
    run:
        sg = swan.read(input.sg)
        die, genes = sg.die_gene_test(obs_col='genotype',
                                      obs_conditions=[wildcards.genotype1,
                                                      wildcards.genotype2])
        die.to_csv(output.out, sep='\t')

rule swan_output_adata:
    input:
        sg = config['data']['sg']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        out = config['data']['adata']
    run:
        save_swan_adata(input.sg,
                        output.out,
                        how=wildcards.feature)


################################################################################
##################################### DEG / DET ################################
################################################################################
rule deg:
    input:
        adata = config['data']['adata']
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = config['data']['de_tsv']
    conda:
        "pydeseq2"
    shell:
        """
            python diff_exp.py \
                   {input.adata} \
                   {wildcards.feature} \
                   genotype \
                   {wildcards.genotype1},{wildcards.genotype2} \
                   {output.out} \
                   {resources.threads}
        """

#         sg = swan.read(input.sg)
#         sg.adata.obs['genotype'] = sg.adata.obs.dataset.str.rsplit('_', n=2, expand=True)[0]
#         die, genes = sg.die_gene_test(obs_col='genotype',
#                                       obs_conditions=[wildcards.genotype1,
#                                                       wildcards.genotype2])
#         die.to_csv(output.out, sep='\t')


################################################################################
########################## BigWig stuff ########################################
################################################################################

def get_strand_flag(wc):
    if wc.strand == 'fwd':
        flag = '--filterRNAstrand reverse'
    elif wc.strand == 'rev':
        flag = '--filterRNAstrand forward'
    return flag

rule sort_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools sort \
            --threads {resources.threads} \
            -O bam {input.bam} > {output.bam}
        """

rule index_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools index -@ {resources.threads} {input.bam}
        """

rule bam_to_bw:
    resources:
        threads = 4,
        mem_gb = 64
    params:
        strand_flag = lambda wc: get_strand_flag(wc)
    conda:
        'deeptools'
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} {params.strand_flag}
        """

def get_input_sr_bam(wc, sr_df):
    return sr_df.loc[sr_df.mouse_id==wc.mouse_id, 'fname'].values[0]

use rule sort_bam as sr_sort_bam with:
    input:
        bam = lambda wc: get_input_sr_bam(wc, sr_df)
    output:
        bam = config['sr']['bam_sorted']

use rule index_bam as sr_index_bam with:
    input:
        bam = config['sr']['bam_sorted']
    output:
        bam = config['sr']['bam_index']

use rule bam_to_bw as sr_bam_to_bw with:
    input:
        bam = config['sr']['bam_sorted'],
        bai = config['sr']['bam_index']
    output:
        bw = config['sr']['bw']

# lr
use rule sort_bam as lr_sort_bam with:
    input:
        bam = config['data']['bam_label_merge']
    output:
        bam = config['data']['bam_label_merge_sorted']

use rule index_bam as lr_index_bam with:
    input:
        bam = config['data']['bam_label_merge_sorted']
    output:
        bam = config['data']['bam_label_merge_index']

use rule bam_to_bw as lr_bam_to_bw with:
    input:
        bam = config['data']['bam_label_merge_sorted'],
        bai = config['data']['bam_label_merge_index']
    output:
        bw = config['data']['bw']


################################################################################
############################ Debugging #########################################
################################################################################

rule subset_bam:
    input:
        bam = config['data']['bam_label_merge_sorted'],
        bai = config['data']['bam_label_merge_index']
    resources:
        threads = 4,
        mem_gb = 16
    params:
        region = 'chr11:102881298-102886826'
    output:
        bam = config['data']['bam_subset']
    shell:
        """
        module load samtools
        samtools view -h {input.bam} {params.region} > {output.bam}
        """

use rule sort_bam as lr_subset_sort_bam with:
    input:
        bam = config['data']['bam_subset']
    output:
        bam = config['data']['bam_subset_sorted']

use rule index_bam as lr_subset_index_bam with:
    input:
        bam = config['data']['bam_subset_sorted']
    output:
        bam = config['data']['bam_subset_index']
