rule talon_label:
    resources:
        threads = 1,
        mem_gb = 32
    shell:
        """
        talon_label_reads \
            --f {input.sam} \
            --g {input.fa} \
            --tmpDir {params.opref} \
            --ar {params.frac_a_range} \
            --deleteTmp \
            --o {params.opref}
        """
