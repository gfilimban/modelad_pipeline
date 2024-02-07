# scratch

# 240203


```python
import cerberus
ref_ics = 'analysis/240202_grant/cerberus/agg/ad003_5xFAD-HEMI_F_4_months_HC_1_ic.tsv'
ics = '/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update_2/data/cerberus/ad003/5xFAD-WT_F_4_months_HC_ic.tsv'
refs = [False, False]
sources = ['cerberus', 'ad003_5xFAD-WT_F_4_months_HC']
output = 'analysis/240202_grant/cerberus/agg/ad003_5xFAD-WT_F_4_months_HC_2_ic.tsv'
cerberus.agg_ics([ref_ics, ics],
                  refs,
                  sources,
                  output)

```

```python
# import cerberus
import pyranges as pr
ref_ends = 'ref/ca_tes.bed'
df = pr.read_bed(ref_ends)
print(df.stranded)

df = df.df
df = pr.PyRanges(df)
print(df.stranded)
# ends='/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update_2/data/cerberus/ad003/5xFAD-HEMI_F_4_months_HC_tes.bed'
# refs = [False for i in range(2)]
# add_ends = [True for i in range(2)]
# cerberus.agg_ends([ref_ends, ends],
#                   add_ends,
#                   refs,
#                   ['cerberus', 's2'],
#                   'tes',
#                   20,
#                   'test_agg.bed')

```

# 240202
```bash
# conda create -n modelad_snakemake python==3.7.5 snakemake
conda install -n base -c conda-forge mamba

mamba create -c conda-forge -c bioconda -n modelad_snakemake_2 snakemake==7.32 python==3.9 pandas pytables==3.8.0

git clone git@github.com:fairliereese/pyranges.git
cd pyranges
pip install -e .
cd ../

pip install sorted-nearest==0.0.37

# conda install pytables==3.8.0

git clone git@github.com:mortazavilab/cerberus.git
cd cerberus
pip install -e .
cd ../

pip install swan-vis

conda install pyfaidx

conda install pandarallel

conda install pysam==0.22.0

conda install pandas
pip install numpy==1.23.5







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
```
