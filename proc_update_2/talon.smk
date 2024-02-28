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

rule talon_init:
	resources:
		mem_gb = 32,
		threads = 16
	shell:
		"""
		talon_initialize_database \
    		--f {input.gtf} \
    		--g {params.genome_ver} \
    		--a {params.annot_ver} \
    		--l {params.min_transcript_len} \
    		--idprefix TALON \
    		--5p {params.max_5_dist} \
    		--3p {params.max_3_dist} \
    		--o {params.opref}
		"""

rule talon:
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
            -v 1
        tmp_dir={params.opref}_temp/
        rm -r ${{tmp_dir}}
        """

rule talon_annot:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_fetch_reads \
            --db {input.db} \
            --build {params.genome_ver} \
            --o {params.opref}
		"""

rule talon_abundance:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_abundance \
            --db {input.db} \
            -a {params.annot_ver} \
            -b {params.genome_ver} \
            --o {params.opref}
		"""

rule talon_filter:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_filter_transcripts \
            --db {input.db} \
            -a {params.annot_ver} \
            --maxFracA {params.max_frac_a}\
            --minCount {params.min_count} \
            --minDatasets {params.min_datasets} \
            --filter_known \
            --o {output.pass_list}
		"""

rule talon_filtered_abundance:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_abundance \
            --db {input.db} \
            -a {params.annot_ver} \
            -b {params.genome_ver} \
            --whitelist {input.pass_list} \
            --o {params.opref}
		"""

rule talon_gtf:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_create_GTF \
            --db {input.db} \
            -a {params.annot_ver} \
            -b {params.genome_ver} \
            --whitelist {input.pass_list} \
            --o {params.opref}
		"""
