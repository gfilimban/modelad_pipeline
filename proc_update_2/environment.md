## Setting up your processing environment

### Analysis pipeline environment

I recommend using [Conda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) to manage your environments. For the analysis portion, the main pipeline is run in the `modelad_snakemake_analysis` environment, and the differential expression tests are run in the `modelad_snakemake_pydeseq2` environment.

To make the `modelad_snakemake_analysis` environment, run the following commands:
```bash
conda install -n base -c conda-forge mamba

mamba create -c conda-forge -c bioconda -n modelad_snakemake_analysis snakemake==7.32 python==3.9 pandas pytables==3.8.0
conda activate modelad_snakemake_analysis

git clone git@github.com:fairliereese/pyranges.git
cd pyranges
pip install -e .
cd ../

pip install sorted-nearest==0.0.37

git clone git@github.com:mortazavilab/cerberus.git
cd cerberus
pip install -e .
cd ../

git clone git@github.com:mortazavilab/swan_vis.git
cd swan_vis
pip install -e .
cd ../

conda install pyfaidx

conda install pandarallel

conda install pysam==0.22.0

conda install pandas
pip install numpy==1.23.5
```

To make the `modelad_snakemake_pydeseq2` environment, run the following commands. *It is important that you use the same name here or the pipeline will not work!*
```bash
mamba create -c conda-forge -c bioconda -n modelad_snakemake_pydeseq2  python==3.9
conda activate modelad_snakemake_pydeseq2
pip install pydeseq2
```



<!-- The YAML files for each are included in this repository. Run the code below to create the environments.
```bash
conda env create -f snakemake.yml
conda env create -f igvtools.yml
conda env create -f deeptools.yml
conda env create -f pydeseq2.yml
conda env create -f snakemake_vis.yml
``` -->

<!-- ```bash
conda activate snakemake
conda env export > snakemake.yml

conda activate igvtools
conda env export > igvtools.yml

conda activate deeptools
conda env export > deeptools.yml

conda activate pydeseq2
conda env export > pydeseq2.yml

conda activate snakemake_vis
conda env export > snakemake_vis.yml
``` -->

<!-- Additionally, you'll need to install the following, in your "snakemake" environment.

* [Cython](https://cython.readthedocs.io/en/latest/src/quickstart/install.html)
```bash
pip install Cython
```
* [samtools](https://anaconda.org/bioconda/samtools)
```bash
conda install bioconda::samtools
```
* [bedtools](https://anaconda.org/bioconda/bedtools)
```bash
conda install bioconda::bedtools
```
* [TranscriptClean](https://github.com/mortazavilab/TranscriptClean)
```bash
git clone git@github.com:mortazavilab/TranscriptClean.git
cd TranscriptClean
pip install -e .
```
* [TALON](https://github.com/mortazavilab/TALON)
```bash
git clone git@github.com:mortazavilab/TALON.git
cd TALON
pip install -e .
```
* [GENCODE UTR Fix](https://github.com/MuhammedHasan/gencode_utr_fix)
```bash
pip install cython
pip install -e git+https://github.com/MuhammedHasan/gencode_utr_fix.git#egg=gencode_utr_fix
```
* [Fairlie's version of PyRanges](git@github.com:fairliereese/pyranges.git)
```bash
git clone git@github.com:fairliereese/pyranges.git
cd pyranges
pip install -e .
``` -->
