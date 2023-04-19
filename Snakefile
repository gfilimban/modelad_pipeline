import pandas as pd
import os
import sys

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from utils import *

configfile: 'config.yml'

df = pd.read_csv('230418_modelad_config.tsv', sep='\t')
batches = df.batch.tolist()
datasets = df.dataset.tolist()
samples = df['sample'].tolist()
platforms = df.platform.tolist()

def get_df_col(wc, df, col):
    val = df.loc[df.dataset==wc.dataset, col].values[0]
    return val

# # config formatting errors
# if len(df.source.unique()) < len(df.source.tolist()):
#     raise ValueError('Sources must have unique names')


rule all:
   #  expand(config['data']['map_stats'],
   #         zip,
   #         batch=batches,
   #         dataset=datasets),
   # expand(config['data']['tc_stats'],
   #        zip,
   #        batch=batches,
   #        dataset=datasets),
   input:
       expand(config['data']['sam_label'],
                      zip,
                      batch=batches,
                      dataset=datasets),
       # config['data']['talon_config']

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
        fa = config['ref']['fa_gz']
    output:
        out = config['ref']['fa']

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
        python {params.tc_path}accessory_scripts/get_SJs_from_gtf.py \
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
        samtools stats {input.alignment} | \
           grep ^SN | \
           cut -f 2- | \
           grep -e 'reads mapped' -e 'reads unmapped' -e 'average length' -e 'maximum length' > {output.stats}
           sed -i '/reads mapped and paired/d' > {output.stats}
        """

################################################################################
################################ Mapping #######################################
################################################################################
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
                --junc-bed {input.sjs}
     		    {input.ref_fa} {input.fastq} > {output.sam} 2> {output.log}"""

rule sam_to_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools sort \
            --threads {params.threads} \
            -O bam {input.sam} > {output.bam}
        samtools index -$ {params.threads} {output.bam}
        """

use rule map as map_reads with:
    input:
        fastq = lambda wc: get_df_col(wc, df, 'fname'),
        ref_fa = config['ref']['fa'],
        sjs = config['ref']['sjs']
    output:
        sam = config['data']['sam'],
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
        sam_rev = config['data']['sam_rev']
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
        sam = config['data']['sam_clean']

use rule alignment_stats as tc_stats with:
    input:
        alignment = config['data']['sam_clean']
    output:
        stats = config['data']['tc_stats']

################################################################################
################################# TALON ########################################
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
        sam = config['data']['sam_label']
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

rule talon_config:
    input:
        files = expand(config['data']['sam_label'],
                       zip,
                       batch=batches,
                       dataset=datasets)
    resources:
        threads = 1,
        mem_gb = 1
    params:
        df = df
    output:
        config = config['data']['talon_config']
    run:
        config = params.df[['dataset', 'sample', 'platform']].copy(deep=True)
        config['fname'] = input.files
        config.to_csv(output.config, header=None, sep=',')
