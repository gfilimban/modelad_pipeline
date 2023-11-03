rule igv_sort_gtf:
    conda:
        "igvtools"
    shell:
        """
        igvtools sort {input.gtf} {output.gtf}
        """

rule igv_index_gtf:
    conda:
        "igvtools"
    shell:
        """
        igvtools index {input.gtf}
        """
