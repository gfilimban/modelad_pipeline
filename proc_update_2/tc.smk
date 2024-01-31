rule get_annot_sjs:
    input:
        gtf = config['ref']['gtf'],
        fa = config['ref']['fa']
    resources:
        threads = 4,
        mem_gb = 16
    params:
        path = config['tc']['path'],
        min_intron_size = config['tc']['min_intron_size']
    output:
        sjs = config['ref']['sjs']
    shell:
        """
        python {params.path}accessory_scripts/get_SJs_from_gtf.py \
             --f {input.gtf} \
             --g {input.fa} \
             --minIntronSize {params.min_intron_size} \
             --o {output.sjs}
        """

rule tc:
    resources:
        mem_gb = 256,
        threads = 16
    shell:
        """
        transcriptclean \
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
