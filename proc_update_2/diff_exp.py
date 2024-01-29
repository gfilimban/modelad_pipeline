from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import scanpy as sc
import sys
import numpy as np

def run_deseq2(fname,
               how,
               obs_col,
               obs_conditions,
               ofile,
               threads):
    """
    Parameters:
        fname (str): Name of input AnnData
        obs_col (str): String column name of condition in obs table
            to test on
        obs_conditions (list of str, len 2): Set of conditions from `obs_col`
            to compare
        how (str): {'gene', 'iso'}
        ofile (str): Name of output file
        threads (int): Number of threads to run on
    """

    # wrong number of conditions
    if len(obs_conditions) != 2:
        raise ValueError(f'{obs_conditions} not valid. Please provide exactly two conditions to compare')

    adata = sc.read(fname)

    # # todo - rm this
    # adata.obs['genotype'] = adata.obs['dataset'].str.rsplit('_', expand=True)[0]

    # column doesn't exist
    if obs_col not in adata.obs.columns:
        raise ValueError(f'{obs_col} not found in adata')
    # remove novel genes
    if how == 'gene':
        adata = adata[:, adata.var.loc[adata.var.index.str.contains('ENSMUS')].index].copy()

    # remove entries that don't belong to the groups we want to test
    adata = adata[adata.obs.loc[adata.obs[obs_col].isin(obs_conditions)].index,
                  :].copy()

    # remove unexpressed stuff
    sc.pp.filter_genes(adata, min_counts=1)

    # densify matrix
    adata.X = np.array(adata.X.todense())

    # run test
    dds = DeseqDataSet(adata=adata,
                   design_factors=obs_col,
                   refit_cooks=True)
    dds.deseq2()
    stat_res = DeseqStats(dds)
    stat_res.summary()

    # save output
    df = stat_res.results_df
    df.to_csv(ofile, sep='\t')

def main():
    fname = sys.argv[1]
    how = sys.argv[2]
    obs_col = sys.argv[3]
    obs_conditions = sys.argv[4].split(',')
    # import pdb
    # pdb.set_trace()
    ofile = sys.argv[5]
    threads = int(sys.argv[6])

    run_deseq2(fname,
                   how,
                   obs_col,
                   obs_conditions,
                   ofile,
                   threads)

if __name__ == '__main__':
    main()
