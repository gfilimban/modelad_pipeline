## Setting up your processing environment

I recommend using [Conda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) to manage your environments.

* [TranscriptClean](https://github.com/mortazavilab/TranscriptClean)
* [TALON](https://github.com/mortazavilab/TALON)
* [GENCODE UTR Fix](https://github.com/MuhammedHasan/gencode_utr_fix)
* [igvtools](https://anaconda.org/bioconda/igvtools)
* [deeptools](https://deeptools.readthedocs.io/en/develop/content/installation.html)
* TODO add yaml config files for enviroemnt
* Get familiar with using `tmux` to stash your sessions while performing preprocessing
* Always use `srun` on HPC3 before running Snakemake or you'll get yelled at


To run this pipeline, I use several Conda environments. The YAML files for each are included in this repository. Run the code below to create the environments.
```bash
conda env create -f snakemake.yml
conda env create -f igvtools.yml
conda env create -f deeptools.yml
conda env create -f pydeseq2.yml
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
``` -->
