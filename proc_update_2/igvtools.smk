rule igv_sort_gtf:
    resources:
        mem_gb = 16,
        threads = 4
    conda:
        "igvtools"
    shell:
        """
        igvtools sort {input.gtf} {output.gtf}
        """

rule igv_index_gtf:
    resources:
        mem_gb = 16,
        threads = 4
    conda:
        "igvtools"
    shell:
        """
        igvtools index {input.gtf}
        """
