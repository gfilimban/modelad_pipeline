# utilities related to defining snakemake rules
from snakemake.io import expand
import pandas as pd
import warnings
import os

def process_meta(meta_fname):
    """
    Process metadat information. Remove spaces and forward
    slashes with underscores. Verify that mouse ids aren't
    duplicated.
    """
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
                      auto_dedupe=True,
                      include_humanized=True):

    """
    Parameters:
        fname (str): Path to config file fname. One line per input fastq.
        meta_fname (str): Path to file with metadata information.
        datasets_per_run (int): Number of datasets to process in each TALON run
        auto_dedupe (bool): Automatically deduplicate duplicate fastqs that result from
            successive Porechop rounds
        include_humanized (bool): Include models with humanized loci, which need
            some preprocessing / different treatment

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

    # source for cerberus should be study + sample
    df['source'] = df['study']+'_'+df['sample']

    # # dataset should be sample + bio rep + flow cel
    # df['dataset'] = df['talon_dataset']+'_'+df['flowcell'].astype(str)

    # get and verify humanized status
    assert len(df.loc[(df.humanized==True)&~(df.genotype.str.contains('h'))]) == 0
    temp = df.loc[(df.humanized==False)&(df.genotype.str.contains('h'))].copy(deep=True)
    if len(temp.index) >= 1:
        genotypes = temp.genotype.unique().tolist()
        warnings.warn(f'Config found non-humanized mouse w/ genotypes {genotypes}, is this expected?')

    if not include_humanized:
        df = df.loc[df.humanized==False].copy(deep=True)

    # assign a cerberus run to each "sample" (study+genotype+sex+age+tissue)
    # but first sort on study and sample such that they will always be ordered in the same way
    # this should freeze our results
    gb_cols = ['study', 'genotype', 'sex', 'age', 'tissue']
    df = df.sort_values(by=gb_cols, ascending=True)
    temp = df.copy(deep=True)
    temp = temp[gb_cols].groupby(gb_cols).count().reset_index()
    temp['cerberus_run'] = [i+1 for i in range(len(temp.index))]
    df = df.merge(temp, how='left', on=gb_cols)

    df['flowcell'] = df.flowcell.astype(str)
    df['biorep_num'] = df.biorep_num.astype(str)
    df['cerberus_run'] = df.cerberus_run.astype(str)

    return df


def subset_df_on_wcs(wc, df):
    """
    Return a copy of the input metadata df limited to the wildcards
    """
    temp = df.copy(deep=True)
    for key, item in wc.items():

        # for entries that we don't have in the df, such as LAPA end mode
        if key not in df.columns:
            continue
        # if key == 'cerberus_run':
        #     import pdb; pdb.set_trace()

        # if we're given a list of possibilities
        if type(item) == list:
            temp = temp.loc[temp[key].isin(item)]
        else:
            temp = temp.loc[temp[key] == item]
    return temp

def get_df_col(wc, df, col):
    """
    From the metadata dataframe df, get the entries that satisfy
    the wildcards requirements and return the corresponding value
    from col. Ensure that this is always a 1:1 relationship, otherwise
    throw an error.
    """
    cols = [col] + [key for key, item in wc.items() if key in df.columns]

    temp = subset_df_on_wcs(wc, df)
    temp = temp[cols].drop_duplicates()

    if len(temp.index) != 1:
        msg = 'Issues getting data from DF with wildcards'
        for key, item in wc.items():
            msg+=f'\n{key}: {item}'
        raise ValueError(msg)

    val = temp[col].tolist()[0]
    return val


def get_cfg_entries(wc, df, cfg_entry, return_df=False):
    """
    Expand a config entry based on the wildcards and the
    values in df that satisfy these wildcards

    Parameters:
        return_df (bool): Return DataFrame with 'file' column
            as opposed to list of files. Default: False
    """
    temp = subset_df_on_wcs(wc, df)

    # cols = ['study', 'genotype', 'sex',
    #         'age', 'tissue', 'biorep_num',
    #         'flowcell']

    study = temp.study.tolist()
    genotype = temp.genotype.tolist()
    sex = temp.sex.tolist()
    age = temp.age.tolist()
    tissue = temp.tissue.tolist()
    biorep_num = temp.biorep_num.tolist()
    flowcell = temp.flowcell.tolist()
    cerberus_run = temp.cerberus_run.tolist()

    files = expand(cfg_entry,
                   zip,
                   study=study,
                   genotype=genotype,
                   sex=sex,
                   age=age,
                   tissue=tissue,
                   biorep_num=biorep_num,
                   flowcell=flowcell,
                   cerberus_run=cerberus_run,
                   allow_missing=True)

    temp['file'] = files

    # make sure we only take unique ones
    # cols.append('file')
    # temp = temp[cols]
    temp = temp.drop_duplicates(subset='file', keep='first')
    files = temp['file'].tolist()

    if return_df:
        return temp
    else:
        return files

def get_lapa_settings(wc, df, cfg_entry, kind):
    """
    Get the command name or file output name
    given whether we're running LAPA in tss or tes mode

    Parameters:
        lapa_ends (str): Config entry
        kind (str): {'temp_file', 'lapa_cmd'}
    """

    # get the output file from the lapa run
    lapa_ends = get_cfg_entries(wc, df, cfg_entry, return_df=False)

    # now expand for end mode
    lapa_ends = expand(lapa_ends,
                       zip=True,
                       end_mode=wc['end_mode'])
    assert len(list(set(lapa_ends))) == 1
    lapa_ends = lapa_ends[0]

    if kind == 'temp_file':
        temp = os.path.dirname(lapa_ends)+'/'
        if wc['end_mode'] == 'tes':
            temp += 'polyA_clusters.bed'
        elif wc['end_mode'] == 'tss':
            temp += 'tss_clusters.bed'
        return temp
    elif kind == 'lapa_cmd':
        if wc['end_mode'] == 'tes':
            return 'lapa'
        elif wc['end_mode'] == 'tss':
            return 'lapa_tss'

def get_first_cerb_entry(wc, df, cfg_entry):
    """
    Get the first config entry run for Cerberus.
    """
    temp_wc = {'cerberus_run': str(1)}
    file = get_cfg_entries(temp_wc, df, cfg_entry)
    assert len(file) == 1
    file = file[0]

    # add in the end mode if we have it
    if 'end_mode' in wc.keys():
        file = expand(file,
                      zip,
                      end_mode=wc['end_mode'])
        assert len(file) == 1
        file = file[0]

    return file

def get_prev_cerb_entry(wc, df, cfg_entry):
    """
    Get the previous config entry run for Cerberus. Ensure that
    only one file meets these criteria.
    """
    prev_wc = {'cerberus_run': str(int(wc['cerberus_run'])-1)}
    file = get_cfg_entries(prev_wc, df, cfg_entry)
    assert len(file) == 1
    file = file[0]

    # add in the end mode if we have it
    if 'end_mode' in wc.keys():
        file = expand(file,
                      zip,
                      end_mode=wc['end_mode'])
        assert len(file) == 1
        file = file[0]

    return file
