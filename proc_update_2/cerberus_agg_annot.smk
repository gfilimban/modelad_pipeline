################################################################################
######################### Cerberus aggregation #################################
################################################################################

rule cerb_agg_ends:
    resources:
        threads = 4,
        mem_gb = 64
    run:
        # only aggregating 2 things at a time
        refs = [params.ref for i in range(2)]
        add_ends = [params.add_ends for i in range(2)]
        cerberus.agg_ends([input.ref_ends, input.ends],
                          add_ends,
                          refs,
                          params.sources,
                          wildcards.end_mode,
                          params.slack,
                          output.ends)

################################################################################
######################### Cerberus annotation ##################################
################################################################################

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

rule cerberus_agg_ics:
    resources:
        mem_gb = 32,
        threads = 2
    run:
        # only aggregating 2 things at a time
        refs = [params.ref for i in range(2)]
        cerberus.agg_ics([input.ref_ics, input.ics],
                          refs,
                          params.sources,
                          output.ics)

rule cerb_write_ref:
  resources:
      threads = 4,
      mem_gb = 64
  run:
      cerberus.write_reference(input.tss,
                               input.tes,
                               input.ic,
                               output.h5)

rule cerb_annot:
   resources:
       mem_gb = 64,
       threads = 16
   run:
       cerberus.annotate_transcriptome(input.gtf,
                                       input.h5,
                                       params.source,
                                       params.gene_source,
                                       output.h5)

rule cerb_gtf_ids:
    resources:
       mem_gb = 64,
       threads = 16
    run:
       cerberus.replace_gtf_ids(input.h5,
                                input.gtf,
                                params.source,
                                params.update_ends,
                                params.agg,
                                output.gtf)

rule cerb_ab_ids:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.replace_ab_ids(input.ab,
                                input.h5,
                                params.source,
                                params.agg,
                                output.ab)
