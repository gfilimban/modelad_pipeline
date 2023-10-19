# utilities related to defining snakemake rules
from snakemake.io import expand

def subset_df_on_wcs(wc, df):
    """
    Return a copy of the input metadata df limited to the wildcards
    """
    temp = df.copy(deep=True)
    for key, item in wc.items():
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
    cols = [col] + [key for key, item in wc.items()]

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

    files = expand(cfg_entry,
                   zip,
                   study=study,
                   genotype=genotype,
                   sex=sex,
                   age=age,
                   tissue=tissue,
                   biorep_num=biorep_num,
                   flowcell=flowcell)

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
