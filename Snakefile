import pandas as pd
import os
import sys
import swan_vis as swan
import itertools

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from utils import *

# settings we can change each time it's run
configfile: 'config.yml'
config_tsv = '230516_config.tsv'
# config_tsv = '230429_config.tsv'
meta_tsv = 'mouse_metadata.tsv'
datasets_per_run = 4 # number of datasets per talon run
auto_dedupe = True # deduplicate runs w/ same stem but different chop numbers

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
studies = temp.study.tolist()
max_talon_runs = temp.talon_run_num.tolist()

wildcard_constraints:
    genotype1= '|'.join([re.escape(x) for x in genotypes]),
    genotype2= '|'.join([re.escape(x) for x in genotypes])


end_modes = ['tss', 'tes']

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

rule all:
    input:
        expand(config['data']['map_stats'],
           zip,
           batch=batches,
           dataset=datasets,
           flowcell=flowcells),
        expand(config['data']['tc_stats'],
          zip,
          batch=batches,
          dataset=datasets,
          flowcell=flowcells),
        expand(config['data']['talon_db'],
          batch=batches,
          study=studies,
          talon_run=max_talon_runs),
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


        # expand(config['data']['lapa_ab'], batch=batches)

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
        bam = temporary(config['data']['bam_label_sorted'])

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

def get_talon_run_files(wc, batches, df, config_entry):
    temp = df.loc[df.talon_run_num == int(wc.talon_run)]
    datasets = temp.dataset.tolist()
    files = expand(config_entry,
                   zip,
                   batch=batches,
                   dataset=datasets)
    return files

rule talon_config:
    input:
        # files = expand(config['data']['sam_label'],
        #                zip,
        #                batch=batches,
        #                dataset=datasets)
        files = lambda wc:get_talon_run_files(wc,
                                              batches,
                                              dataset_df,
                                              config['data']['bam_label_merge'])
    resources:
        threads = 1,
        mem_gb = 1
    params:
        df = dataset_df
    output:
        config = config['data']['talon_config']
    run:
        config = params.df[['dataset', 'sample', 'platform', 'talon_run_num']].copy(deep=True)
        config = config.loc[config.talon_run_num==int(wildcards.talon_run)]
        config.drop('talon_run_num', axis=1, inplace=True)
        config['fname'] = input.files
        config.to_csv(output.config, header=None, sep=',', index=False)

rule talon_init:
	input:
		ref_gtf = config['ref']['gtf']
	output:
		db = config['ref']['talon_db']
	params:
		talon_opref = config['ref']['talon_db'].rsplit('_talon', maxsplit=1)[0],
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
        cp {input.ref} {input.ref}_back
        talon \
            --f {input.config} \
            --db {input.ref} \
            --build {params.genome} \
            --tmpDir {params.opref}_temp/ \
            --threads {resources.threads} \
            --o {params.opref}
        mv {input.ref} {params.opref}_talon.db
        mv {input.ref}_back {input.ref}
        """

use rule talon as first_talon with:
    input:
        ref = config['ref']['talon_db'],
        config = expand(config['data']['talon_config'],
                        batch=batches,
                        talon_run=1)[0]
    params:
        genome = 'mm10',
        opref = expand(config['data']['talon_db'],
                        batch=batches,
                        talon_run=1)[0].rsplit('_talon', maxsplit=1)[0],
    output:
        db = expand(config['data']['talon_db'],
                        batch=batches,
                        talon_run=1)[0],
        annot = expand(config['data']['read_annot'],
                        batch=batches,
                        talon_run=1)[0]

use rule talon as seq_talon with:
    input:
        ref = lambda wc: expand(config['data']['talon_db'],
                        batch=wc.batch,
                        talon_run=int(wc.talon_run)-1)[0],
        config = config['data']['talon_config']
    params:
        genome = 'mm10',
        opref = config['data']['talon_db'].rsplit('_talon', maxsplit=1)[0]
    output:
        db = config['data']['talon_db'],
        annot = config['data']['read_annot']

rule talon_unfilt_ab:
    input:
        db = expand(config['data']['talon_db'],
                        batch=batches,
                        talon_run=max_talon_run)[0]
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
        db = expand(config['data']['talon_db'],
                        batch=batches,
                        talon_run=max_talon_run)[0]
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
        db = expand(config['data']['talon_db'],
                        batch=batches,
                        talon_run=max_talon_run)[0],
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
        db = expand(config['data']['talon_db'],
                        batch=batches,
                        talon_run=max_talon_run)[0],
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
################################ Swan ##########################################
################################################################################
def make_sg(input, params, wildcards):

    # initialize
    sg = swan.SwanGraph()
    sg.add_annotation(input.annot)
    sg.add_transcriptome(input.gtf, include_isms=True)
    sg.add_abundance(input.filt_ab)
    sg.add_abundance(input.ab, how='gene')

    # add metadata
    adatas = sg.get_adatas()
    for adata in adatas:
        adata.obs['genotype'] = adata.obs.dataset.str.rsplit('_', n=2, expand=True)[0]

    sg.save_graph(params.prefix)

rule make_sg:
    input:
        annot = config['ref']['gtf'],
        gtf = config['data']['filt_gtf'],
        filt_ab = config['data']['filt_ab'],
        ab = config['data']['ab']
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
################################# LAPA #########################################
################################################################################
def get_lapa_settings(wc, lapa_ends, kind):
    lapa_ends = expand(lapa_ends, zip, batch=wc.batch, end_mode=wc.end_mode)[0]
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

rule lapa_config:
    input:
        # todo - this will likely become config['data']['bam_label_merge']
        # files = expand(config['data']['sam_label'],
        #                zip,
        #                batch=batches,
        #                dataset=datasets)
    resources:
        threads = 1,
        mem_gb = 1
    params:
        df = df
    output:
        config = config['data']['lapa_config']
    run:
        config = params.df[['sample', 'dataset']].copy(deep=True)
        config['fname'] = input.files
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
        annot = config['data']['read_annot'],
        tss = expand(config['data']['lapa_ends'], end_mode='tss', batch=batches)[0],
        tes = expand(config['data']['lapa_ends'], end_mode='tes', batch=batches)[0]
    resources:
        threads = 1,
        mem_gb = 64
    params:
        tss_dir = expand(config['data']['lapa_ends'], end_mode='tss', batch=batches)[0].rsplit('/', maxsplit=1)[0]+'/',
        tes_dir = expand(config['data']['lapa_ends'], end_mode='tes', batch=batches)[0].rsplit('/', maxsplit=1)[0]+'/'
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
        annot = config['data']['read_annot'],
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
        mem_gb = 4
    params:
        t_nov = ['Known', 'NIC', 'NNC', 'ISM_rescue'],
        g_nov = ['Known'],
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


rule filt_lapa_gtf:


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
