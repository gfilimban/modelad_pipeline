import cerberus

end_modes = ['tss', 'tes']

################################################################################
####################### Cerberus references ####################################
################################################################################

use rule dl as dl_ca with:
    params:
        link = config['ref']['cerberus']['link']
    output:
        out = config['ref']['cerberus']['ca']

rule ca_to_tsv:
    input:
        ca = config['ref']['cerberus']['ca']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        opref = config['ref']['cerberus']['ics'].split('_ic.tsv')[0]
    output:
        expand(config['ref']['cerberus']['ends'],
               zip,
               end_mode=end_modes),
        config['ref']['cerberus']['ics']
    run:
        cerberus.write_h5_to_tsv(input.ca,
                                 params.opref)

################################################################################
#################### Get triplet features from GTF #############################
################################################################################
rule cerb_gtf_to_bed:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_bed(input.gtf,
                            wildcards.end_mode,
                            output.ends,
                            dist=params.dist,
                            slack=params.slack)

rule cerb_gtf_to_ics:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_ics(input.gtf,
                            output.ics)

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
