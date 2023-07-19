import yaml
from snakemake.io import expand
import swan_vis as swan
import pandas as pd
import numpy as np

def process_meta(meta_fname):
    meta = pd.read_csv(meta_fname, sep='\t')
    meta['genotype'] = meta.genotype.str.replace(' ', '_')
    meta['genotype'] = meta.genotype.str.replace('/', '_')

    meta['age'] = meta.age.str.replace(' ', '_')

    dupe_ids = meta.loc[meta.mouse_id.duplicated(), 'mouse_id'].tolist()
    if len(dupe_ids) > 0:
        raise ValueError(f'Found duplicated mouse ids {dupe_ids} in mouse metadata')

    return meta

def parse_config_file(fname,
                      meta_fname,
                      datasets_per_run,
                      auto_dedupe=True):

    """
    Parameters:
        fname (str): Path to config file fname. One line per input fastq.
        meta_fname (str): Path to file with metadata information.
        datasets_per_run (int): Number of datasets to process in each TALON run
        auto_dedupe (bool): Automatically deduplicate duplicate fastqs that result from
            successive Porechop rounds

    Returns:
        df (pandas DataFrame): DF w/ pipeline information; one line per fastq
        dataset_df (pandas DataFrame): DF w/ dataset information; one line per mouse
    """

    df = pd.read_csv(fname, sep='\t')

    ############ Basename + fname
    df['basename'] = df.fname.str.rsplit('/', n=1, expand=True)[1]
    df['path'] = df.fname.str.rsplit('/', n=1, expand=True)[0]

    ############ Dataset + flowcell df

    # get flowcell
    exp = '.*\/[\w-]+_(\d+)(?:_t\d+)?\.fastq(?:.gz)?'
    df['flowcell'] = df.fname.str.extract(exp)

    # check to make sure the same file stem isn't there more than once
    # (can happen if different flow cells needed different amounts of chopping)
    # df['file_stem'] = df.basename.str.rsplit('_', n=1, expand=True)[0]
    exp = '.*\/([\w-]+_\d+)(?:_t\d+)?\.fastq(?:.gz)?'
    df['file_stem'] = df.fname.str.extract(exp)
    df['chop_num'] = df.basename.str.rsplit('.fastq', expand=True)[0].str.rsplit('_t', expand=True)[1].astype(float)
    if df.file_stem.duplicated().any():
        dupe_stems = df.loc[df.file_stem.duplicated(keep=False), 'basename'].tolist()
        if not auto_dedupe:
            raise ValueError(f'Files {dupe_stems} seem to be duplicated. Check config file.')
        else:
            print(f'Files {dupe_stems} seem to be duplicated. Automatically removing lower chop numbers')
            df = df.sort_values(by='chop_num', ascending=False)
            df = df.drop_duplicates(subset='file_stem', keep='first')

    # extract the sample name
    temp = df.basename.str.split('_', expand=True)[[0,1]]#.str.join('_')
    df['sample_temp'] = temp[0]+'_'+temp[1]

    # extract the mouse id
    df['mouse_id'] = df['sample_temp'].str.split('_', expand=True)[1]

    # extract the "study" name
    exp = '^(ad[0-9]+)'
    df['study'] = df.basename.str.extract(exp)


    # merge in metadata
    meta = process_meta(meta_fname)
    df['mouse_id'] = df['mouse_id'].astype('int')
    df = df.merge(meta, how='left', on='mouse_id')

    # get tech rep numbers -- each mouse has multiple reps
    # and are therefore technical reps
    df['flowcell'] = df.sort_values(['genotype', 'mouse_id'],
                                ascending=[True, True])\
                                .groupby(['mouse_id']) \
                                .cumcount() + 1

    # sample should be the genotype + age + sex + tissue
    df['sample'] = df.genotype+'_'+ \
                   df.sex+'_'+ \
                   df.age+'_'+ \
                   df.tissue

    # get biorep numbers -- each mouse_id is a different mouse
    # and therefore a different biorep
    temp = df[['sample', 'mouse_id']].drop_duplicates()
    temp.reset_index(inplace=True, drop=True)
    temp['biorep_num'] = temp.sort_values(['sample', 'mouse_id'],
                                ascending=[True, True])\
                                .groupby(['sample']) \
                                .cumcount()+1
    df = df.merge(temp, how='left',
                  on=['sample', 'mouse_id'])

    # talon dataset should be sample + bio rep
    df['dataset'] = df['sample']+'_'+df['biorep_num'].astype(str)

    # # dataset should be sample + bio rep + flow cel
    # df['dataset'] = df['talon_dataset']+'_'+df['flowcell'].astype(str)

    ############ TALON dataset df

    # create a dataset-level df that will represent the aggregate
    cols = ['sample', 'mouse_id', 'genotype', 'sex', 'study', \
            'age', 'tissue', 'biorep_num', 'dataset', 'platform']
    dataset_df = df[cols].drop_duplicates()

    # get the talon run number these will go into
    dataset_df = dataset_df.sort_values(by='study', ascending=False)
    dataset_df = dataset_df.reset_index(drop=True)
    dataset_df['talon_run_num'] = np.nan
    for ind, entry in dataset_df.iterrows():
        curr_study = entry.study
        # first entry
        if ind == 0:
            run_num = 0
            study_ind = 0
            prev_study = curr_study

        # when we find a new study
        if curr_study != prev_study:
            run_num = 0
            study_ind = 0

        # when we hit the max. # datasets / run
        if study_ind % datasets_per_run == 0:
            run_num += 1

        # actual assignment of number
        dataset_df.loc[ind, 'talon_run_num'] = run_num

        # keep track of last study that was assigned a number
        prev_study = curr_study

        # inc study ind
        study_ind += 1

    # clean up some typing
    dataset_df['talon_run_num'] = dataset_df.talon_run_num.astype(int)
    df['flowcell'] = df.flowcell.astype(str)

    return df, dataset_df


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
new_meta = 'swan_metadata.tsv'
meta2 = pd.read_csv(meta, sep='\t')
meta2.to_csv(new_meta, sep='\t', index=False) # og ad005
meta2 = pd.read_csv(ad003_meta, sep='\t')
meta2.to_csv(new_meta, sep='\t', index=False, mode='a') # ad003

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
