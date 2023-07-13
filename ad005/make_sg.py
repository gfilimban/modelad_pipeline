import yaml
from snakemake.io import expand

config_file = '../config.yml'
with open(config_file) as f:
    config = yaml.safe_load(f)

study = 'ad005'
batch = '230516'

annot = expand(config['data']['ca_ref_gtf'],
               zip,
               study=study,
               batch=batch)
gtf = expand(config['data']['ca_gtf'],
               zip,
               study=study,
               batch=batch)
filt_ab = expand(config['data']['ca_ab'],
               zip,
               study=study,
               batch=batch)
ab = expand(config['data']['ab'],
               zip,
               study=study,
               batch=batch)
meta = config['data']['swan_meta']
prefix = config['data']['sg'].replace('.p', '')
sg = expand(config['data']['sg'],
               zip,
               study=study,
               batch=batch)

ad003_gtf = expand(config['data']['ca_gtf'],
               zip,
               study='ad003',
               batch=batch)
ad003_filt_ab =expand(config['data']['ca_ab'],
               zip,
               study='ad003',
               batch=batch)
ad003_ab =expand(config['data']['ab'],
               zip,
               study='ad003',
               batch=batch)


# initialize
sg = swan.SwanGraph()
sg.add_annotation(annot)
sg.add_transcriptome(gtf, include_isms=True)
sg.add_transcriptome(ad003_gtf, include_isms=True)
sg.add_abundance(filt_ab)
sg.add_abundance(ad003_filt_ab)
sg.add_abundance(ab, how='gene')
sg.add_abundance(ad003_ab, how='gene')
sg.save_graph(prefix)

sg.add_metadata(meta)
sg.save_graph(prefix)
