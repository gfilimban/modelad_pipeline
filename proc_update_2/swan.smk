import swan_vis as swan

rule make_swan_metadata:
    params:
        p_df = p_df
    resources:
        mem_gb = 1,
        threads = 1
    output:
        meta = config['analysis']['swan']['meta']
    run:
        temp_meta = get_cfg_entries_analysis(wildcards, params.p_df,
                                    config['merge']['sort_bam'],
                                    return_df=True)
        cols = ['dataset', 'mouse_id', 'study', 'genotype',
                'sex', 'genotype_sex', 'genotype_alias', 'age', 'tissue', 'biorep_num']
        temp_meta = temp_meta[cols]
        temp_meta.to_csv(output.meta, sep='\t', index=False)

rule make_swan_graph:
    input:
        gtf = lambda wc:get_cfg_entries_analysis(wc, p_df, config['analysis']['cerberus']['gtf']),
        ab = lambda wc:get_cfg_entries_analysis(wc, p_df, config['analysis']['cerberus']['ab']),
        annot = config['analysis']['ref']['cerberus']['gtf'],
        meta = config['analysis']['swan']['meta']
    resources:
        mem_gb = 64,
        threads = 4
    params:
        prefix = config['analysis']['swan']['swan_graph'].replace('.p', '')
    output:
        sg = config['analysis']['swan']['swan_graph']
    run:
        sg = swan.SwanGraph()
        sg.add_annotation(input.annot)
        for g in input.gtf:
            sg.add_transcriptome(g, include_isms=True)

        for ab in input.ab:
            sg.add_abundance(ab)

        for ab in input.ab:
            sg.add_abundance(ab, how='gene')

        sg.add_metadata(input.meta)
        sg.save_graph(params.prefix)

rule swan_die:
    input:
        sg = config['analysis']['swan']['swan_graph']
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = config['analysis']['swan']['du']
    run:
        sg = swan.read(input.sg)
        die, genes = sg.die_gene_test(obs_col=wildcards.obs_col,
                                      obs_conditions=[wildcards.obs_cond1,
                                                      wildcards.obs_cond2],
                                      kind=wildcards.feat)
        die.to_csv(output.out, sep='\t')

def save_swan_adata(swan_file,
                    ofile,
                    how='iso'):
    """
    Save anndata obj from Swan.

    Input:
        swan_file (str): Input SwanGraph file
        ofile (str): Output AnnData file
        how (str): {'iso', 'tss', 'tes', 'ic', 'edge', 'loc', 'gene'}
    """
    sg = swan.read(swan_file)
    if how == 'gene':
        adata = sg.gene_adata
    elif how == 'iso':
        adata = sg.adata
    else:
        raise ValueError("You haven't implemented this yet.")
    adata.write(ofile)

rule swan_output_g_adata:
    input:
        sg = config['analysis']['swan']['swan_graph']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        out = temporary(config['analysis']['swan']['g_adata'])
    run:
        save_swan_adata(input.sg,
                        output.out,
                        how='gene')

rule swan_output_t_adata:
    input:
        sg = config['analysis']['swan']['swan_graph']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        out = temporary(config['analysis']['swan']['t_adata'])
    run:
        save_swan_adata(input.sg,
                        output.out,
                        how='iso')

################################################################################
##################################### DEG / DET ################################
################################################################################
rule deg:
    input:
        adata = config['analysis']['swan']['g_adata']
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = config['analysis']['swan']['deg']
    conda:
        "pydeseq2"
    shell:
        """
            python diff_exp.py \
                   {input.adata} \
                   gene \
                   {wildcards.obs_col} \
                   {wildcards.obs_cond1},{wildcards.obs_cond2} \
                   {output.out} \
                   {resources.threads}
        """

rule det:
    input:
        adata = config['analysis']['swan']['t_adata']
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = config['analysis']['swan']['det']
    conda:
        "pydeseq2"
    shell:
        """
            python diff_exp.py \
                   {input.adata} \
                   iso \
                   {wildcards.obs_col} \
                   {wildcards.obs_cond1},{wildcards.obs_cond2} \
                   {output.out} \
                   {resources.threads}
        """

rule all_swan:
    input:
        expand(config['analysis']['swan']['swan_graph'],
               analysis=p_df.analysis.dropna().unique().tolist()),
        get_de_cfg_entries(p_df, config['analysis']['swan']['du'], how='du'),
        get_de_cfg_entries(p_df, config['analysis']['swan']['deg'], how='de'),
        get_de_cfg_entries(p_df, config['analysis']['swan']['det'], how='de')
