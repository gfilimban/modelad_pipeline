rule sam_to_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools view -hSb {input.sam} > {output.bam}
        """

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

rule merge_alignment:
    resources:
        threads = 32,
        mem_gb = 64
    shell:
        """
        module load samtools
        samtools merge -o {output.bam} {input.files}
        """
