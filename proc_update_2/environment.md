## Setting up your processing environment

I recommend using [Conda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) to manage your environments. I use several Conda environments. The YAML files for each are included in this repository. Run the code below to create the environments.
```bash
conda env create -f snakemake.yml
conda env create -f igvtools.yml
conda env create -f deeptools.yml
conda env create -f pydeseq2.yml
conda env create -f snakemake_vis.yml
```

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

Additionally, you'll need to install the following, in your "snakemake" environment.
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
```
