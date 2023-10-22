from utils import *

rule lapa_link:
    resources:
        threads = 1,
        mem_gb = 256
    shell:
        """
        lapa_link_tss_to_tes \
            --alignment {input.annot} \
            --lapa_dir {params.tes_dir} \
            --lapa_tss_dir {params.tss_dir} \
            --output {output.links}
        """

rule lapa_call_ends:
    resources:
        threads = 4,
        mem_gb = 32
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

rule lapa_correct_talon:
    resources:
        threads = 1,
        mem_gb = 128
    shell:
        """
        lapa_correct_talon \
                --links {input.links} \
                --read_annot {input.annot} \
                --gtf_input {input.gtf} \
                --gtf_output {output.gtf} \
                --abundance_input {input.filt_ab} \
                --abundance_output {output.filt_ab} \
                --keep_unsupported
        """

rule lapa_filt:
    resources:
        threads = 1,
        mem_gb = 32
    run:
        filter_lapa(input.filt_ab,
                    input.gtf,
                    params.t_novs,
                    params.g_novs,
                    params.filt_spikes,
                    output.filt_list)

rule lapa_filt_ab:
    resources:
        threads = 4,
        mem_gb = 32
    run:
        df = filt_lapa_ab(input.ab,
                          input.filt_list)
        df.to_csv(output.ab, sep='\t', index=False)

rule lapa_filt_gtf:
    resources:
        threads = 4,
        mem_gb = 32
    run:
        gtf = filt_lapa_gtf(input.gtf,
                            input.filt_list)
        gtf.to_gtf(output.gtf)
