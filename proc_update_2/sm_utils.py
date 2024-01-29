# utilities related to defining snakemake rules
from snakemake.io import expand
import pandas as pd
import warnings
import os
import itertools

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
                      p_meta_fname,
                      auto_dedupe=True):

    """
    Parameters:
        fname (str): Path to config file fname. One line per input fastq.
        meta_fname (str): Path to file with metadata information.
        p_meta_fname (str): Path to pseudochromosome metadata information
        datasets_per_run (int): Number of datasets to process in each TALON run
        auto_dedupe (bool): Automatically deduplicate duplicate fastqs that result from
            successive Porechop rounds
        include_pseudochrom (bool): Include models with pseudochrom loci, which need
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
    # assert len(df.loc[(df.pseudochrom_needed==True)&~(df.genotype.str.contains('h'))]) == 0
    temp = df.loc[(df.pseudochrom_needed==False)&(df.genotype.str.contains('h'))].copy(deep=True)
    # if len(temp.index) >= 1:
    #     genotypes = temp.genotype.unique().tolist()
    #     warnings.warn(f'Config found non-pseudochrom mouse w/ genotypes {genotypes}, is this expected?')

    # if not include_pseudochrom:
    #     df = df.loc[df.pseudochrom_needed==False].copy(deep=True)

    # get pseudochromosome / reference : genotype df
    # else:
    def format_pseudochrom_cols(df, col):
        """
        Format hgene, mgene, and pseudochromosome names columns
        to either replace NaNs with "dummy" and and to string
        split entries with more than one
        """
        inds = df.loc[(df[col].isnull())].index
        df.loc[inds, col] = 'dummy'
        df[col] = df[col].str.split(',')
        df[col] = df.apply(lambda x: tuple(sorted(x[col])), axis=1)
        return df

    # for c in ['pseudochrom', 'human_gene', 'mouse_gene']:
    for c in ['pseudochrom']:
        df = format_pseudochrom_cols(df, c)

    # make sure the correspondance between
    # genotype:pseudochromosomes is 1:1
    temp = df.loc[df.pseudochrom_needed==True].copy(deep=True)
    temp = temp[['pseudochrom', 'genotype']].drop_duplicates()
    dupe_genotypes = temp.loc[temp.genotype.duplicated()].genotype.unique().tolist()
    if len(dupe_genotypes) > 1:
        raise ValueError(f'Found genotype(s) {dupe_genotypes} w/ multiple pseudochromosome settings')

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

    # add in columns for comparisons
    df['genotype_sex'] = df['genotype']+'_'+df['sex']

    # get a table that matches genotype + pseudochrom + human gene + mouse gene
    temp = df.explode('pseudochrom')
    p_meta = pd.read_csv(p_meta_fname, sep='\t')
    p_meta.fillna('dummy', inplace=True)
    p_df = temp.merge(p_meta, on='pseudochrom')

    return df, p_df

def parse_analysis_config_file(fname,
                      meta_fname,
                      p_meta_fname,
                      an_meta_fname,
                      auto_dedupe=True):

    """
    Parameters:
        fname (str): Path to config file fname. One line per input fastq.
        meta_fname (str): Path to file with metadata information.
        p_meta_fname (str): Path to pseudochromosome metadata information
        an_meta_fname (str): Path to analysis metadata information
        datasets_per_run (int): Number of datasets to process in each TALON run
        auto_dedupe (bool): Automatically deduplicate duplicate fastqs that result from
            successive Porechop rounds
        include_pseudochrom (bool): Include models with pseudochrom loci, which need
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
    # assert len(df.loc[(df.pseudochrom_needed==True)&~(df.genotype.str.contains('h'))]) == 0
    temp = df.loc[(df.pseudochrom_needed==False)&(df.genotype.str.contains('h'))].copy(deep=True)
    if len(temp.index) >= 1:
        genotypes = temp.genotype.unique().tolist()
        warnings.warn(f'Config found non-pseudochrom mouse w/ genotypes {genotypes}, is this expected?')

    # if not include_pseudochrom:
    #     df = df.loc[df.pseudochrom_needed==False].copy(deep=True)

    # get pseudochromosome / reference : genotype df
    # else:
    def format_pseudochrom_cols(df, col):
        """
        Format hgene, mgene, and pseudochromosome names columns
        to either replace NaNs with "dummy" and and to string
        split entries with more than one
        """
        inds = df.loc[(df[col].isnull())].index
        df.loc[inds, col] = 'dummy'
        df[col] = df[col].str.split(',')
        df[col] = df.apply(lambda x: tuple(sorted(x[col])), axis=1)
        return df

    # for c in ['pseudochrom', 'human_gene', 'mouse_gene']:
    for c in ['pseudochrom']:
        df = format_pseudochrom_cols(df, c)

    # make sure the correspondance between
    # genotype:pseudochromosomes is 1:1
    temp = df.loc[df.pseudochrom_needed==True].copy(deep=True)
    temp = temp[['pseudochrom', 'genotype']].drop_duplicates()
    dupe_genotypes = temp.loc[temp.genotype.duplicated()].genotype.unique().tolist()
    if len(dupe_genotypes) > 1:
        raise ValueError(f'Found genotype(s) {dupe_genotypes} w/ multiple pseudochromosome settings')

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

    # add in columns for comparisons
    df['genotype_sex'] = df['genotype']+'_'+df['sex']

    # get a table that matches genotype + pseudochrom + human gene + mouse gene
    temp = df.explode('pseudochrom')
    p_meta = pd.read_csv(p_meta_fname, sep='\t')
    p_meta.fillna('dummy', inplace=True)
    p_df = temp.merge(p_meta, on='pseudochrom')

    # add in analysis stuff
    an_df = pd.read_csv(an_meta_fname, sep='\t')
    p_df = p_df.merge(an_df, how='left',
                  on=['genotype', 'study'])

    return df, p_df



def subset_df_on_wcs(wc, df):
    """
    Return a copy of the input metadata df limited to the wildcards
    """
    temp = df.copy(deep=True)
    for key, item in wc.items():

        # for entries that we don't have in the df, such as LAPA end mode
        if key not in df.columns:
            continue

        # if we're given a list of possibilities
        if type(item) == list:
            temp = temp.loc[temp[key].isin(item)]
        else:
            temp = temp.loc[temp[key] == item]
        # if len(temp.index) == 0:
        #     import pdb; pdb.set_trace()
    return temp

def get_df_col(wc, df, col, allow_multiple=False):
    """
    From the metadata dataframe df, get the entries that satisfy
    the wildcards requirements and return the corresponding value
    from col. Ensure that this is always a 1:1 relationship, otherwise
    throw an error.

    Parameters:
        allow_multiple (bool): Whether to allow multiple rows per
            wcs given. Default = False
    """
    cols = [col] + [key for key, item in wc.items() if key in df.columns]
    cols = list(set(cols))

    temp = subset_df_on_wcs(wc, df)
    temp = temp[cols].drop_duplicates()

    if not allow_multiple:
        if len(temp.index) != 1:
            msg = 'Issues getting data from DF with wildcards'
            for key, item in wc.items():
                msg+=f'\n{key}: {item}'
            raise ValueError(msg)
        val = temp[col].tolist()[0]
    else:
        val = temp[col].tolist()

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

    genotype = temp.genotype.tolist()
    study = temp.study.tolist()
    sex = temp.sex.tolist()
    age = temp.age.tolist()
    tissue = temp.tissue.tolist()
    biorep_num = temp.biorep_num.tolist()
    flowcell = temp.flowcell.tolist()
    cerberus_run = temp.cerberus_run.tolist()
    pseudochrom = temp.pseudochrom.tolist()
    # try:
    mouse_gene = temp.mouse_gene.tolist()
    human_gene = temp.human_gene.tolist()
    analysis = temp.analysis.tolist()
    # except:
    #     import pdb; pdb.set_trace()

    # # pseudochrom stuff needs to be treated differently
    # pseudochrom = temp.pseudochrom.tolist()
    # todo need to figure out if this is a correct assertion
    # i'm a lil more confident about it now
    # assert len(set(pseudochrom)) == 1
    # pseudochrom = list(pseudochrom[0])
    #
    # files = expand(cfg_entry,
    #                zip,
    #                study=study,
    #                genotype=genotype,
    #                sex=sex,
    #                age=age,
    #                tissue=tissue,
    #                biorep_num=biorep_num,
    #                flowcell=flowcell,
    #                cerberus_run=cerberus_run,
    #                pseudochrom=pseudochrom,
    #                allow_missing=True)

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
                   pseudochrom=pseudochrom,
                   human_gene=human_gene,
                   mouse_gene=mouse_gene,
                   analysis=analysis,
                   allow_missing=True)
    # import pdb; pdb.set_trace()

    # # if we're working with multiple files for one
    # # entry (ie pseudochroms), turn that into a tuple
    # if len(files) > 1 and len(temp.index) == 1:
    #     files = [tuple(files)]
    #     untuple = True
    # else:
    #     untuple = False

    # try:
    temp['file'] = files
    # except:
    #     import pdb; pdb.set_trace()


    # make sure we only take unique ones
    temp = temp.drop_duplicates(subset='file', keep='first')
    files = temp['file'].tolist()

    if return_df:
        return temp
    else:
    #     if untuple:
    #         assert len(files) == 1
    #         files = list(files[0])
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
    first_cerb_run = df.cerberus_run.min(axis=0)
    temp_wc = {'cerberus_run': str(first_cerb_run)}
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

def get_prev_cerb_entry(wc, df, cfg_entry, config):
    """
    Get the previous config entry run for Cerberus. Ensure that
    only one file meets these criteria.
    """
    prev_run = str(int(wc['cerberus_run'])-1)

    # for first entry, we should be using the reference
    # set of ends / ics
    if prev_run == '0':
        if 'end_mode' in wc.keys():
            file = config['ref']['cerberus']['ends']
        else:
            file = config['ref']['cerberus']['ics']
    else:
        prev_wc = {'cerberus_run': prev_run}
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

def get_de_cfg_entries(p_df, cfg_entry, how):
    """
    Get file names needed as output for DE or DU
    tests within analysis objects for
        - all pairwise genotype sets
        - all pairwise genotype sets by sex

    Parameters:
        how (str): {'du', 'de'}
    """

    if how == 'du':
        feats = ['tss', 'tes', 'ic', 'iso']
    else:
        feats = []

    files = []

    for a in p_df.analysis.unique().tolist():
        # print()
        # print(a)
        wc = {'analysis': a}
        temp = subset_df_on_wcs(wc, p_df)
        obs_col = 'genotype'
        conds = temp[obs_col].unique().tolist()
        # genotypes += ['ghost_cookie']

        combos = [c for c in itertools.combinations(conds, 2)]
        obs_cond1 = [c[0] for c in combos]
        obs_cond2 = [c[1] for c in combos]
        # print(genotype1)
        # print(genotype2)

        files += expand(expand(cfg_entry,
          zip,
          obs_cond1=obs_cond1,
          obs_cond2=obs_cond2,
          allow_missing=True),
          obs_col=obs_col,
          feat=feats,
          analysis=a)

        # now get genotype comparisons for each sex
        for s in temp.sex.unique():
            wc['sex'] = s
            temp2 = subset_df_on_wcs(wc, temp)
            obs_col = 'genotype_sex'
            conds = temp2[obs_col].unique().tolist()

            combos = [c for c in itertools.combinations(conds, 2)]
            obs_cond1 = [c[0] for c in combos]
            obs_cond2 = [c[1] for c in combos]

            files += expand(expand(cfg_entry,
              zip,
              obs_cond1=obs_cond1,
              obs_cond2=obs_cond2,
              allow_missing=True),
              obs_col=obs_col,
              feat=feats,
              analysis=a)

    return files
