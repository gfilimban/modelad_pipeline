# scratch

# 240202
# conda create -n modelad_snakemake python==3.7.5 snakemake
mamba create -c conda-forge -c bioconda -n modelad_snakemake snakemake

git clone git@github.com:fairliereese/pyranges.git
cd pyranges
pip install -e .
cd ../

pip install sorted-nearest==0.0.37

git clone git@github.com:mortazavilab/cerberus.git
cd cerberus
pip install -e .
cd ../




# snakemake env
conda create -n modelad_snakemake python==3.7

conda activate modelad_snakemake





git clone git@github.com:mortazavilab/TALON.git
cd TALON
pip install -e .
cd ../

git clone git@github.com:fairliereese/pyranges.git
cd pyranges
pip install -e .
cd ../

conda install -c conda-forge bioconda::snakemake

conda install bioconda::samtools

conda install bioconda::bedtools

pip install cython

pip install -e git+https://github.com/MuhammedHasan/gencode_utr_fix.git#egg=gencode_utr_fix

pip install sorted-nearest==0.0.33

git clone git@github.com:mortazavilab/cerberus.git
cd cerberus
pip install -e .
cd ../


# transcriptclean environments
git clone git@github.com:mortazavilab/TranscriptClean.git
cd TranscriptClean
pip install -e .
cd ../



# pydeseq2 thing
# pydeseq2



# trying to run lapa separately
lapa_tss \
  --alignment data/lapa/ad004/hTREM2KI-WT_M_8_months_HC/config.csv \
  --fasta ref/genome/hTREM2KI-WT/merged.fa \
  --annotation ref/annot_utr.gtf \
  --chrom_sizes ref/genome/hTREM2KI-WT/chrom_sizes.txt \
  --output_dir data/lapa/ad004/hTREM2KI-WT_M_8_months_HC/tss/
