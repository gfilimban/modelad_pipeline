# make a fasta file concatenating the original reference and the
# pseudochromosomes. also, just symlink the original reference
# if we have a dummy pseudochrom
rule mkref_cat:
    resources:
        threads = 1,
        mem_gb = 4
    run:
        # dummy chr -- just symlink original
        # fasta in the directory for this genotype
        pseudochroms =  get_df_col(wildcards,
                                   params.p_df,
                                   'pseudochrom',
                                    allow_multiple=True)
        if pseudochroms == ['dummy']:
            os.symlink(os.path.abspath(input.ref), output.out)
        # otherwise cat everything together
        else:
            infiles = [input.ref]+input.files
            with open(output.out, 'w') as outfile:
                for fname in infiles:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)

rule tc_human:
    resources:
        mem_gb = 256,
        threads = 16
    shell:
        """
        if [ {wildcards.human_gene} == "dummy" ]
        then
            touch {output.sam}
            touch {output.fa}
            touch {output.sam_clean_log}
            touch {output.sam_clean_te_log}
        else
            python {params.path}TranscriptClean.py \
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
        fi
        """

rule tc_mouse:
    resources:
        mem_gb = 256,
        threads = 16
    shell:
        """
        if [ {wildcards.mouse_gene} == "dummy" ]
        then
            touch {output.sam}
            touch {output.fa}
            touch {output.sam_clean_log}
            touch {output.sam_clean_te_log}
        else
            python {params.path}TranscriptClean.py \
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
        fi
        """

rule sam_to_bam_mouse:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        if [ {wildcards.mouse_gene} == "dummy" ]
        then
            touch {output.bam}
        else
            module load samtools
            samtools view -hSb {input.sam} > {output.bam}
        fi
        """

rule sam_to_bam_human:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        if [ {wildcards.human_gene} == "dummy" ]
        then
            touch {output.bam}
        else
            module load samtools
            samtools view -hSb {input.sam} > {output.bam}
        fi
        """

rule sort_bam_mouse:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        if [ {wildcards.mouse_gene} == "dummy" ]
        then
            touch {output.bam}
        else
            module load samtools
            samtools sort \
                --threads {resources.threads} \
                -O bam {input.bam} > {output.bam}
        fi
        """

rule sort_bam_human:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        if [ {wildcards.human_gene} == "dummy" ]
        then
            touch {output.bam}
        else
            module load samtools
            samtools sort \
                --threads {resources.threads} \
                -O bam {input.bam} > {output.bam}
        fi
        """

rule index_bam_mouse:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        if [ {wildcards.mouse_gene} == "dummy" ]
        then
            touch {output.ind}
        else
            module load samtools
            samtools index -@ {resources.threads} {input.bam}
        fi
        """

rule index_bam_human:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        if [ {wildcards.human_gene} == "dummy" ]
        then
            touch {output.ind}
        else
            module load samtools
            samtools index -@ {resources.threads} {input.bam}
        fi
        """





use rule mkref_cat as mkref_genome with:
    input:
        ref = config['ref']['fa'],
        files = lambda wc: get_cfg_entries(wc,
                                           p_df,
                                           config['ref']['pseudochrom']['fa'])
    params:
        p_df = p_df
    output:
        out = config['ref']['pseudochrom']['fa_merge']

rule mkref_human_gene_fq:
   input:
       fa = config['human_ref']['t_fa']
   resources:
       mem_gb = 16,
       threads = 4
   output:
       fastq = config['ref']['pseudochrom']['human_gene']['fq']
   run:
       get_gene_t_fastq(input.fa,
                        wildcards.human_gene,
                        output.fastq)

use rule map as map_reads_hgene with:
    input:
        fastq = config['ref']['pseudochrom']['human_gene']['fq'],
        ref_fa = config['ref']['pseudochrom']['fa'],
        sjs = config['ref']['sjs']
    output:
        sam = temporary(config['ref']['pseudochrom']['human_gene']['sam']),
        log = config['ref']['pseudochrom']['human_gene']['log']

use rule tc_human as tc_sam_hgene with:
    input:
        sam = rules.map_reads_hgene.output.sam,
        fa = config['ref']['pseudochrom']['fa']
    params:
        path = config['tc']['path'],
        min_intron_size = config['tc']['min_intron_size'],
        opref = config['ref']['pseudochrom']['human_gene']['tc_sam'].rsplit('_clean.sam', maxsplit=1)[0]
    output:
        sam = temporary(config['ref']['pseudochrom']['human_gene']['tc_sam']),
        fa = temporary(config['ref']['pseudochrom']['human_gene']['tc_fa']),
        sam_clean_log = temporary(config['ref']['pseudochrom']['human_gene']['tc_log']),
        sam_clean_te_log = temporary(config['ref']['pseudochrom']['human_gene']['te_log'])

use rule sam_to_bam_human as bam_from_sam_hgene with:
    input:
        sam = rules.tc_sam_hgene.output.sam
    output:
        bam = temporary(config['ref']['pseudochrom']['human_gene']['bam'])

use rule sort_bam_human as bam_sort_hgene with:
    input:
        bam = rules.bam_from_sam_hgene.output.bam
    output:
        bam = config['ref']['pseudochrom']['human_gene']['sort_bam']

use rule index_bam_human as bam_ind_hgene with:
    input:
        bam = rules.bam_sort_hgene.output.bam
    output:
        ind = config['ref']['pseudochrom']['human_gene']['ind_bam']

rule mkref_mouse_gene_fq:
   input:
       fa = config['ref']['t_fa']
   resources:
       mem_gb = 16,
       threads = 4
   output:
       fastq = config['ref']['pseudochrom']['gene']['fq']
   run:
       get_gene_t_fastq(input.fa,
                        wildcards.mouse_gene,
                        output.fastq)

use rule map as map_reads_mgene with:
  input:
      fastq = config['ref']['pseudochrom']['gene']['fq'],
      ref_fa = config['ref']['pseudochrom']['fa'],
      sjs = config['ref']['sjs']
  output:
      sam = config['ref']['pseudochrom']['gene']['sam'],
      log = config['ref']['pseudochrom']['gene']['log']

use rule tc_mouse as tc_sam_mgene with:
  input:
      sam = rules.map_reads_mgene.output.sam,
      fa = config['ref']['pseudochrom']['fa']
  params:
      path = config['tc']['path'],
      min_intron_size = config['tc']['min_intron_size'],
      opref = config['ref']['pseudochrom']['gene']['tc_sam'].rsplit('_clean.sam', maxsplit=1)[0]
  output:
      sam = temporary(config['ref']['pseudochrom']['gene']['tc_sam']),
      fa = temporary(config['ref']['pseudochrom']['gene']['tc_fa']),
      sam_clean_log = temporary(config['ref']['pseudochrom']['gene']['tc_log']),
      sam_clean_te_log = temporary(config['ref']['pseudochrom']['gene']['te_log'])

use rule sam_to_bam_mouse as bam_from_sam_mgene with:
  input:
      sam = rules.tc_sam_mgene.output.sam
  output:
      bam = temporary(config['ref']['pseudochrom']['gene']['bam'])

use rule sort_bam_mouse as bam_sort_mgene with:
  input:
      bam = rules.bam_from_sam_mgene.output.bam
  output:
      bam = config['ref']['pseudochrom']['gene']['sort_bam']

use rule index_bam_mouse as bam_ind_mgene with:
  input:
      bam = rules.bam_sort_mgene.output.bam
  output:
      ind = config['ref']['pseudochrom']['gene']['ind_bam']

########################################################
############ For pseudochromosome mappings #############
########################################################
rule talon_pseudochrom:
  resources:
      mem_gb = 256,
      threads = 30
  shell:
      """
      cp {input.db} {output.db}
      talon \
          --f {input.cfg} \
          --db {output.db} \
          --build {params.genome_ver} \
          --tmpDir {params.opref}_temp/ \
          --threads {resources.threads} \
          --create_novel_spliced_genes \
          --o {params.opref} \
          --identity 0.0 \
          --cov 0.0 \
          -v 1
      """

def mk_pseudochrom_mapped_gene_talon_config(wc,
                                          species,
                                          p_df,
                                          ofile):
      if species == 'mouse':
          cfg_entry = config['ref']['pseudochrom']['gene']['sort_bam']
          col = 'mouse_gene'
      elif species == 'human':
          cfg_entry = config['ref']['pseudochrom']['human_gene']['sort_bam']
          col = 'human_gene'
      gene = get_df_col(wc, p_df, col)

      if gene == 'dummy':
          pathlib.Path(ofile).touch()

      else:
          # input files
          bam = get_cfg_entries(wc,
                                p_df,
                                cfg_entry)

          # actual config file makery
          temp = get_cfg_entries(wc,
                                 p_df,
                                 cfg_entry,
                                 return_df=True)
          temp = temp[['sample', 'dataset', 'platform', 'file']]
          temp['sample'] = gene
          temp['dataset'] = gene
          temp.to_csv(ofile, index=False, header=None, sep=',')

rule talon_config_pseudochrom_mouse:
    input:
        bam_ind = rules.bam_ind_mgene.output.ind
    resources:
        threads = 1,
        mem_gb = 2
    params:
        species = 'mouse'
    output:
        cfg = config['ref']['pseudochrom']['gene']['config']
    run:
        mk_pseudochrom_mapped_gene_talon_config(wc,
                                              params.species,
                                              p_df,
                                              output.cfg)

rule talon_config_pseudochrom_human:
  input:
      bam_ind = rules.bam_ind_hgene.output.ind
  resources:
      threads = 1,
      mem_gb = 2
  params:
      species = 'human'
  output:
      cfg = config['ref']['pseudochrom']['human_gene']['config']
  run:
      mk_pseudochrom_mapped_gene_talon_config(wc,
                                            params.species,
                                            p_df,
                                            output.cfg)


rule talon_gtf_pseudochrom:
  resources:
      mem_gb = 16,
      threads = 1
  shell:
      """
      talon_create_GTF \
          --db {input.db} \
          -a {params.annot_ver} \
          -b {params.genome_ver} \
          --observed \
          --o {params.opref}
      """


# use rule talon_pseudochrom as talon_pseudochrom_mouse with:
#     input:
#         db = config['ref']['talon']['db'],
#         config
#     output:
#         db = config['ref']['pseudochrom']['gene']['db']

rule all_pseudochrom:
    input:
        list(set(expand(config['ref']['pseudochrom']['gene']['config'],
               zip,
               pseudochrom=p_df.pseudochrom.tolist(),
               mouse_gene=p_df.mouse_gene.tolist()))),
        list(set(expand(config['ref']['pseudochrom']['human_gene']['config'],
               zip,
               pseudochrom=p_df.pseudochrom.tolist(),
               human_gene=p_df.human_gene.tolist())))
