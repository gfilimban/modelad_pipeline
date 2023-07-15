import yaml
from snakemake.io import expand
import swan_vis as swan

config_file = '../config.yml'
with open(config_file) as f:
    config = yaml.safe_load(f)

study = 'ad005'
batch = '230516'

annot = '../'+expand(config['data']['ca_ref_gtf'],
               zip,
               study=study,
               batch=batch)[0]
gtf = '../'+expand(config['data']['ca_gtf'],
               zip,
               study=study,
               batch=batch)[0]
filt_ab = '../'+expand(config['data']['ca_ab'],
               zip,
               study=study,
               batch=batch)[0]
ab = '../'+expand(config['data']['ab'],
               zip,
               study=study,
               batch=batch)[0]
meta = '../'+expand(config['data']['swan_meta'],
                    zip,
                    study=study,
                    batch=batch)[0]
# prefix = config['data']['sg'].replace('.p', '')
prefix = 'test'
sg = '../'+expand(config['data']['sg'],
               zip,
               study=study,
               batch=batch)[0]

ad003_gtf = '../'+expand(config['data']['ca_gtf'],
               zip,
               study='ad003',
               batch=batch)[0]
ad003_filt_ab ='../'+expand(config['data']['ca_ab'],
               zip,
               study='ad003',
               batch=batch)[0]
ad003_ab ='../'+expand(config['data']['ab'],
               zip,
               study='ad003',
               batch=batch)[0]
ad003_meta = '../'+expand(config['data']['swan_meta'],
                   zip,
                   study='ad003',
                   batch=batch)[0]


# get concatenated metadata file
config_tsv = '230516_config.tsv'
cerb_tsv = 'cerberus.tsv' # gtf_to_bed and agg_ends settings
datasets_per_run = 4 # number of datasets per talon run
auto_dedupe = True # deduplicate runs w/ same stem but different chop numbers
cerb_settings = pd.read_csv(cerb_tsv, sep='\t')
# should be static
meta_tsv = 'mouse_metadata.tsv'
df, dataset_df = parse_config_file(config_tsv,
                       meta_tsv,
                       datasets_per_run=datasets_per_run,
                       auto_dedupe=auto_dedupe)
dataset_df.to_csv(meta, sep='\t', index=False, mode='a')

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
