################################################################################
####################### Cerberus references ####################################
################################################################################

use rule dl as dl_ca with:
    params:
        link = config['ref']['ca_link']
    output:
        out = config['ref']['ca']

rule ca_to_tsv:
    input:
        ca = config['ref']['ca']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        opref = config['ref']['ca_ics'].split('_ic.tsv')[0]
    output:
        expand(config['ref']['ca_ends'],
               zip,
               end_mode=end_modes),
        config['ref']['ca_ics']
    run:
        cerberus.write_h5_to_tsv(input.ca,
                                 params.opref)
