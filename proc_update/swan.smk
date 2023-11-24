import swan_vis as swan

rule make_swan_metadata:
    params:
        p_df = p_df
    resources:
        mem_gb = 1,
        threads = 1
    output:
        meta = confqig['swan']['meta']
    run:
        temp_meta = get_cfg_entries(wildcards, params.p_df,
                                    config['merge']['sort_bam'],
                                    return_df=True)
        cols = ['dataset', 'mouse_id', 'study', 'genotype',
                'sex', 'age', 'tissue', 'biorep_num']
        temp_meta = temp_meta[cols]
        temp_meta.to_csv(output.meta, sep='\t', index=False)

rule make_swan_graph:
    input:
        gtf = lambda wc:get_cfg_entries(wc, p_df, config['cerberus']['gtf']),
        ab = lambda wc:get_cfg_entries(wc, p_df, config['cerberus']['ab']),
        annot = config['ref']['cerberus']['gtf'],
        meta = config['swan']['meta']
    resources:
        mem_gb = 64,
        threads = 4
    params:
        prefix = config['swan']['swan_graph'].replace('.p', '')
    output:
        sg = config['swan']['swan_graph']
    run:
        sg = swan.SwanGraph()
        sg.add_annotation(input.annot)
        for g in input.gtf:
            sg.add_transcriptome(g, include_isms=True)

        for ab in input.ab:
            sg.add_abundance(ab)

        sg.add_metadata(input.meta)
        sg.save_graph(params.prefix)

rule all_swan:
    input:
        expand(config['swan']['swan_graph'],
               analysis=p_df.analysis.dropna().unique().tolist())
