use rule cerb_agg_ends as cerb_agg_ends_lr with:
    input:
        ref_ends = lambda wc: get_prev_cerb_entry(wc, df, config['cerberus']['agg']['ends']),
        ends = config['cerberus']['ends']
    params:
        add_ends = True,
        ref = False,
        slack = lambda wc:config['cerberus']['agg'][wc.end_mode]['slack'],
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ends = config['cerberus']['agg']['ends']


############## old cerberus code

use rule cerberus_agg_ics as cerb_agg_ics_first with:
    input:
        ref_ics = config['ref']['cerberus']['ics'],
        ics = lambda wc:get_first_cerb_entry(wc, df, config['cerberus']['ics'])
    params:
        ref = False,
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ics = config['cerberus']['agg']['ics']

use rule cerberus_agg_ics as cerb_agg_ics_seq with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, df, config['cerberus']['agg']['ics']),
        ics = config['cerberus']['ics']
    params:
        ref = False,
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ics = config['cerberus']['agg']['ics']

use rule cerb_agg_ends as cerb_agg_ends_first with:
    input:
        ref_ends = config['ref']['cerberus']['ends'],
        ends = lambda wc:get_first_cerb_entry(wc, df, config['cerberus']['ends'])
    params:
        add_ends = True,
        ref = False,
        slack = lambda wc:config['cerberus']['agg'][wc.end_mode]['slack'],
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ends = config['cerberus']['agg']['ends']

use rule cerb_agg_ends as cerb_agg_ends_seq with:
    input:
        ref_ends = lambda wc: get_prev_cerb_entry(wc, df, config['cerberus']['agg']['ends']),
        ends = config['cerberus']['ends']
    params:
        add_ends = True,
        ref = False,
        slack = lambda wc:config['cerberus']['agg'][wc.end_mode]['slack'],
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        ends = config['cerberus']['agg']['ends']



######################################
rule cerberus_agg_ics_cfg:
    resources:
        threads = 1,
        mem_gb = 1
    run:
        refs = [params.ref for i in range(2)]
        df = pd.DataFrame()
        df['fname'] = [input.ref_ics, input.ics]
        df['ref'] = refs
        df['sources'] = params.sources
        df.to_csv(output.cfg, sep=',', header=None, index=False)

rule cerberus_agg_ics_cli:
    resources:
        threads = 2,
        mem_gb = 32
    shell:
        """
        cerberus agg_ics \
            --input {input.cfg} \
            -o {output.ics}
        """

use rule cerberus_agg_ics_cfg as cerb_agg_ics_cfg_lr with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['analysis']['cerberus']['agg']['ics'],
                                                  config,
                                                  p_dir),
        ics = p_dir+config['cerberus']['ics']
    params:
        ref = False,
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        cfg = config['analysis']['cerberus']['agg']['ics_cfg']

use rule cerberus_agg_ics_cli as cerb_agg_ics_cfg_cli_lr with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['analysis']['cerberus']['agg']['ics'],
                                                  config,
                                                  p_dir),
        ics = p_dir+config['cerberus']['ics'],
        cfg = config['analysis']['cerberus']['agg']['ics_cfg']
    output:
        ics = config['analysis']['cerberus']['agg']['ics']
