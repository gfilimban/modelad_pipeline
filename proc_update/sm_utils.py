# utilities related to defining snakemake rules

def get_df_col(wc, df, col):
    """
    From the metadata dataframe df, get the entries that satisfy
    the wildcards requirements and return the corresponding value
    from col. Ensure that this is always a 1:1 relationship, otherwise
    throw an error.
    """
    temp = df.copy(deep=True)
    cols = [col]
    for key, item in wc.items():
        temp = temp.loc[temp[key] == item]
        cols.append(key)
    temp = temp[cols].drop_duplicates()
    
    if len(temp.index) != 1:
        msg = 'Issues getting data from DF with wildcards'
        for key, item in wc.items():
            msg+=f'\n{key}: {item}'
        raise ValueError(msg)

    val = temp[col].tolist()[0]
    return val
