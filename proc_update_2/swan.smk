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
                'sex', 'genotype_sex', 'genotype_alias',
                'genotype_alias_int', 'age', 'tissue', 'biorep_num']
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

def plot_du_plot(df, wc, params, ofile):
    from adjustText import adjust_text
    import matplotlib.pylab as plt
    import numpy as np

    df = pd.read_csv(df, sep='\t')

    df['DU'] = False
    df.loc[(df.adj_p_val<=params.adj_p_thresh)&\
           (df.dpi.abs()>=params.dpi_thresh), 'DU'] = True

    # counts of each
    num_du = len(df.loc[df.DU].index)
    num_not_du = len(df.loc[~df.DU].index)

    # add pseudocount to significance and -log10
    df['sig'] = -np.log10(df.adj_p_val+0.01)

    plt.scatter(x=df['dpi'], y=df['sig']., s=1,
                label=f'Not significant (n={num_not_du})')
    du = df.loc[df.DU]
    du.sort_values(['sig'], inplace=True)
    plt.scatter(x=du['dpi'], y=df['sig'], s=3,
                label=f"DU b/w {wc['obs_cond1']} and {wc['obs_cond2']} (n={num_du})", color='b')
    texts = []
    for i in range(min(10, du.shape[0])):
        texts.append(plt.text(x=du.iloc[i]['dpi'],
                              y=du.iloc[i]['sig'],
                              s=du.iloc[i]['gname']))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel("logFC")
    plt.ylabel("-log10(adj p-value+0.01)")
    plt.axvline(0, color="grey", linestyle="--")
    plt.axhline(-np.log10(0.05), color="grey", linestyle="--")
    plt.legend(loc='upper right', fontsize=7)

    plt.savefig(ofile, dpi=500)

rule du_plot:
    input:
        du = config['analysis']['swan']['du']['du']
    params:
        adj_p_thresh = config['analysis']['swan']['du']['adj_p_thresh'],
        dpi_thresh = config['analysis']['swan']['du']['dpi_thresh']
    resources:
        mem_gb = 8,
        threads = 1
    output:
        fname = config['analysis']['swan']['du']['du_plot']
    run:
        plot_du_plot(input.du, wildcards, params, output.fname)

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
def filt_de(sg, de, params, ofile, kind='gene'):
    """
    Format DE table and filter based on thresholds. Add gene names.
    """

    df = pd.read_csv(de, sep='\t')

    # add gene names
    sg = swan.read(sg)
    sg.t_df = sg.t_df.reset_index(drop=True)
    if kind == 'transcript':
        g_df = sg.t_df[['tid', 'tname']].drop_duplicates().reset_index()
        merge_thing = 'tid'
    elif kind == 'gene':
        g_df = sg.t_df[['gid', 'gname']].drop_duplicates().reset_index()
        df['gid_stable'] = cerberus.get_stable_gid(df, 'gid')
        df.drop('gid', axis=1, inplace=True)
        df.rename({'gid_stable':'gid'}, axis=1, inplace=True)
        merge_thing = 'gid'
    df = df.merge(g_df, how='left', on=merge_thing)

    # call things as upregulated or downregulated
    df['DE'] = 'No'
    df.loc[(df.log2FoldChange >= params.l2fc_thresh)&\
           (df.padj <= params.adj_p_thresh), 'DE'] = 'Up'
    df.loc[(df.log2FoldChange <= -1*params.l2fc_thresh)&\
           (df.padj <= params.adj_p_thresh), 'DE'] = 'Down'

    df.to_csv(ofile, sep='\t', index=False)

def plot_v_plot(df, wc, params, ofile, kind='gene'):
    from adjustText import adjust_text
    import matplotlib.pylab as plt
    import numpy as np

    df = pd.read_csv(df, sep='\t')

    if kind == 'gene':
        df['label'] = df.gname
    elif kind == 'transcript':
        df['label'] = df.tname

    if kind == 'gene':
        label = 'gname'
    elif kind == 'transcript':
        label = 'tname'

    df.loc[df.DE == "No", 'label'] = ""
    # Calculate counts
    num_up = df[df.DE == "Up"].shape[0]
    num_down = df[df.DE == "Down"].shape[0]
    num_not_significant = df[df.DE == "No"].shape[0]
    # Plotting
    plt.scatter(x=df['log2FoldChange'], y=df['padj'].apply(lambda x: -np.log10(x)), s=1,
                label=f"Not significant (n={num_not_significant})")
    down = df[df.DE == "Down"]
    down.sort_values(["padj"], inplace=True)
    plt.scatter(x=down['log2FoldChange'], y=down['padj'].apply(lambda x: -np.log10(x)), s=3,
                label=f"Down-regulated in {wc['obs_cond1']} (n={num_down})", color="blue")
    up = df[df.DE == "Up"]
    up.sort_values(["padj"], inplace=True)
    plt.scatter(x=up['log2FoldChange'], y=up['padj'].apply(lambda x: -np.log10(x)), s=3,
                label=f"Up-regulated in {wc['obs_cond1']} (n={num_up})", color="red")
    texts = []
    for i in range(min(10, up.shape[0])):
        texts.append(plt.text(x=up.iloc[i]['log2FoldChange'],
                              y=-np.log10(up.iloc[i]['padj']),
                              s=up.iloc[i][label]))
    for i in range(min(10, down.shape[0])):
            texts.append(plt.text(x=down.iloc[i]['log2FoldChange'],
                              y=-np.log10(up.iloc[i]['padj']),
                              s=down.iloc[i][label]))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel("logFC")
    plt.ylabel("-log10(adj p-value)")
    plt.axvline(params.l2fc_thresh, color="grey", linestyle="--")
    plt.axhline(-np.log10(params.adj_p_thresh), color="grey", linestyle="--")
    # Adjust the legend with a numerical font size
    plt.legend(loc='upper right', fontsize=7)  # Change the font size here

    plt.savefig(ofile, dpi=500)

rule deg:
    input:
        adata = config['analysis']['swan']['g_adata']
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = temporary(config['analysis']['swan']['deg']['deg'])
    conda:
        "modelad_snakemake_pydeseq2"
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

rule deg_fmt:
    input:
        de = config['analysis']['swan']['deg']['deg'],
        sg = config['analysis']['swan']['swan_graph']
    params:
        l2fc_thresh = config['analysis']['swan']['deg']['l2fc_thresh'],
        adj_p_thresh = config['analysis']['swan']['deg']['adj_p_thresh']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        fname = config['analysis']['swan']['deg']['deg_fmt']
    run:
        filt_de(input.sg, input.de, params, output.fname)

rule deg_plot:
    input:
        degs = config['analysis']['swan']['deg']['deg_fmt'],
    resources:
        mem_gb = 64,
        threads = 1
    params
        l2fc_thresh = config['analysis']['swan']['deg']['l2fc_thresh'],
        adj_p_thresh = config['analysis']['swan']['deg']['adj_p_thresh']
    output:
        fname = config['analysis']['swan']['deg']['deg_plot']
    run:
        plot_v_plot(input.degs, wildcards, params, output.fname)

rule det:
    input:
        adata = config['analysis']['swan']['t_adata']
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = temporary(config['analysis']['swan']['det']['det'])
    conda:
        "modelad_snakemake_pydeseq2"
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

rule det_fmt:
    input:
        de = config['analysis']['swan']['det']['det'],
        sg = config['analysis']['swan']['swan_graph']
    params:
        l2fc_thresh = config['analysis']['swan']['det']['l2fc_thresh'],
        adj_p_thresh = config['analysis']['swan']['det']['adj_p_thresh']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        fname = config['analysis']['swan']['det']['det_fmt']
    run:
        filt_de(input.sg, input.de, params, output.fname, kind='transcript')

rule det_plot:
    input:
        dets = config['analysis']['swan']['det']['det_fmt'],
    resources:
        mem_gb = 64,
        threads = 1
    params
        l2fc_thresh = config['analysis']['swan']['det']['l2fc_thresh'],
        adj_p_thresh = config['analysis']['swan']['det']['adj_p_thresh']
    output:
        fname = config['analysis']['swan']['det']['det_plot']
    run:
        plot_v_plot(input.dets, wildcards, params, output.fname, kind='transcript')

rule all_swan:
    input:
        expand(config['analysis']['swan']['swan_graph'],
               analysis=p_df.analysis.dropna().unique().tolist()),
        get_de_cfg_entries(p_df, config['analysis']['swan']['deg']['deg_plot'], how='de'),
        get_de_cfg_entries(p_df, config['analysis']['swan']['det']['det_plot'], how='de'),
        get_de_cfg_entries(p_df, config['analysis']['swan']['du']['du_plot'], how='du'),
