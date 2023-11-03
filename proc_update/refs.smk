import os
# configfile: 'config.yml'

# pseudochrom = config['params']['pseudochrom']

################################################################################
######################## Ref. download / proc ##################################
################################################################################

use rule dl as dl_human_t_fa with:
    params:
        link = lambda wc:config['human_ref']['t_fa_link']
    output:
        out = temporary(config['human_ref']['t_fa_gz'])

use rule gunzip as gz_human_t_fa with:
    input:
        gz = config['human_ref']['t_fa_gz']
    output:
        out = config['human_ref']['t_fa']

use rule dl as dl_t_fa with:
    params:
        link = lambda wc:config['ref']['t_fa_link']
    output:
        out = temporary(config['ref']['t_fa_gz'])

use rule gunzip as gz_t_fa with:
    input:
        gz = config['ref']['t_fa_gz']
    output:
        out = config['ref']['t_fa']

use rule dl as dl_fa with:
    params:
        link = lambda wc:config['ref']['fa_link']
    output:
        out = temporary(config['ref']['fa_gz'])

use rule gunzip as gz_fa with:
    input:
        gz = config['ref']['fa_gz']
    output:
        out = config['ref']['fa']

use rule dl as dl_annot with:
    params:
        link = lambda wc:config['ref']['gtf_link']
    output:
        out = temporary(config['ref']['gtf_gz'])

use rule gunzip as gz_annot with:
    input:
        gz = config['ref']['gtf_gz']
    output:
        out = config['ref']['gtf']

# rule mkref_spike_gtf:
#     input:
#         sirv = config['ref']['spike']['sirv_gtf'],
#         ercc = config['ref']['spike']['ercc_gtf'],
#         gtf = config['ref']['gtf']
#     resources:
#         mem_gb = 16,
#         threads = 2
#     output:
#         all = config['ref']['talon']['gtf']
#     shell:
#         """
#         cat {input.gtf} > {output.all}
#         cat {input.sirv} >> {output.all}
#         cat {input.ercc} >> {output.all}
#         """
#
# rule mkref_spike_fa:
#     input:
#         sirv = config['ref']['spike']['sirv_fa'],
#         ercc = config['ref']['spike']['ercc_fa'],
#         fa = config['ref']['fa']
#     resources:
#         threads = 1,
#         mem_gb = 4
#     output:
#         cat_fa = config['ref']['talon']['fa']
#     shell:
#         """
#         cat {input.fa} >> {output.cat_fa}
#         cat {input.sirv} >> {output.cat_fa}
#         cat {input.ercc} >> {output.cat_fa}
#         """

# make a fasta file concatenating the original reference and the
# pseudochromosomes. also, just symlink the original reference
# if we have a dummy pseudochrom
rule mkref_cat_fastas:
    resources:
        threads = 1,
        mem_gb = 4
    run:
        # dummy chr -- just symlink original
        # fasta in the directory for this genotype
        if wildcards.pseudochrom == 'dummy':
            os.symlink(input.fa, output.fa)
        # otherwise cat everything together
        else:
            infiles = [input.fa]+input.files
            with open(output.fa, 'w') as outfile:
                for fname in infiles:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)

        # """
        # cat {input.fa1} >> {output.fa}
        # cat {input.fa2} >> {output.fa}
        # """

rule mkref_chrom_sizes:
    input:
        fa = config['ref']['pseudochrom']['fa_merge']
    resources:
        threads = 1,
        mem_gb = 8
    output:
        chrom_sizes = config['ref']['pseudochrom']['chrom_sizes']
    shell:
        """
        faidx {input.fa} -i chromsizes > {output.chrom_sizes}
        """

# installable from here: https://github.com/MuhammedHasan/gencode_utr_fix
rule mkref_utr_fix_gtf:
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

################################################################################
######################## Human ref. download ###################################
################################################################################

use rule dl as dl_human_fa with:
    params:
        link = lambda wc:config['human_ref']['fa_link']
    output:
        out = temporary(config['human_ref']['fa_gz'])

use rule gunzip as gz_human_fa with:
    input:
        gz = config['human_ref']['fa_gz']
    output:
        out = config['human_ref']['fa']

use rule dl as dl_human_annot with:
    params:
        link = lambda wc:config['human_ref']['gtf_link']
    output:
        out = temporary(config['human_ref']['gtf_gz'])

use rule gunzip as gz_human_annot with:
    input:
        gz = config['human_ref']['gtf_gz']
    output:
        out = config['human_ref']['gtf']


# ################################################################################
# ########################### Pseudochromosome ###################################
# ################################################################################
#
# use rule mkref_cat_fas as mkref_pseudochrom with:
#     input:
#         fa1 = config['ref']['fa'],
#         fa2 = expand(config['ref']['pseudochrom']['fa'],
#                     pseudochrom=pseudochrom)[0]
#     output:
#         fa = expand(config['ref']['pseudochrom']['fa_merge'],
#                    pseudochrom=pseudochrom)[0]
