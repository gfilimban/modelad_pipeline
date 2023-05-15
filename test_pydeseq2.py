import swan_vis as swan
import pydeseq2
import pdb
import scanpy as sc

swan_file = '/share/crsp/lab/seyedam/freese/mortazavi_lab/bin/modelad_pipeline/data/230429/swan/swan.p'
how='gene'
threads=8
obs_col='genotype'
obs_conditions=['5xFADHEMI', '5xCLU-h2kbKI_HO']

"""
Parameters:
    obs_col (str): String column name of condition in obs table
        to test on
    obs_conditions (list of str, len 2): Set of conditions from `obs_col`
        to compare
    how (str): {'gene', 'iso'}
"""
pdb.set_trace()

# wrong number of condition
if len(obs_conditions) != 2:
    raise ValueError(f'{obs_conditions} not valid. Please provide exactly two conditions to compare')

sg = swan.read(swan_file)
if how == 'gene':
    adata = sg.gene_adata
elif how == 'iso':
    adata = sg.adata

# column doesn't exist
if obs_col not in adata.obs.columns:
    raise ValueError(f'{obs_col} not found in adata')

# remove entries that don't belong to the groups we want to test
adata = adata[adata.obs.loc[adata.obs[obs_col].isin(obs_conditions)].index,
              :].copy(deep=True)

# remove unexpressed stuff
adata = sc.pp.filter_genes(adata, min_counts=1)

# create deseq obj
dds = pydeseq2.dds.DeseqDataSet(adata=adata,
                                design_factors=obs_col,
                                n_cpus=threads,
                                refit_cooks=True)

# run test
dds.deseq2()
stat_res = pydeseq2.ds.DeseqStats(dds,
                                  n_cpus=threads)
# results = stat_res.s
