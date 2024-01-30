# MODEL-AD Long-read RNA-seq pipeline


## Preprocessing
Preparing your processing environment:
* TranscriptClean
* gencode_utr_fix
* igvtools
* deeptools
* TODO add yaml config files for enviroemnt
* TODO mention qi and tmux


* When you get new reads, add the paths to FASTQ files to `config.tsv`.
* When you get new mice, add the relevant data
* `Snakefile` will run the preprocessing pipeline through LAPA.

Run the preprocessing pipeline with the following command (replace `{your_email}` with your email you want updates at!):
```bash
snakemake \
-s Snakefile \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A \
  model-ad_lab \
  --partition=highmem \
  --mem={resources.mem_gb}GB \
  -c {resources.threads} \
  --mail-user={your_email} \
  --mail-type=START,END,FAIL \
  --time=72:00:00" \
  -n
```

Wait for this command to run and make sure the steps that it plans to run are reasonable. After it finishes, run the same command without the `-n` option.

## Downstream analysis

For the majority of you, you will only be running this part of the pipeline. This part of the pipeline will run some basic analyses on the long-read data for you including differential gene and transcript expression tests and differential isoform usage tests.

* First, make sure you *clone (`git clone`) and / or update (`git pull`) this repo* so that you have access to all the relevant metadata.
* When you want to run a new analysis, edit or create a new `analysis_config.tsv`. Check the `mouse_metadata.tsv` file for the relevant genotypes and studies to include in your table.
* `Snakefile_analysis.smk` will run the preprocessing pipeline through LAPA.

Run the preprocessing pipeline with the following command (replace `{your_email}` with your email you want updates at!):
```bash
snakemake \
-s Snakefile_analysis.smk \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A \
  model-ad_lab \
  --partition=highmem \
  --mem={resources.mem_gb}GB \
  -c {resources.threads} \
  --mail-user={your_email} \
  --mail-type=START,END,FAIL \
  --time=72:00:00" \
  -n
```

<!-- ```bash
snakemake \
-s Snakefile_analysis.smk \
-j 200 \
--latency-wait 120 \
--use-conda \
  -n

snakemake \
-s Snakefile_analysis.smk \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A \
  model-ad_lab \
  --partition=highmem \
  --mem={resources.mem_gb}GB \
  -c {resources.threads} \
  --mail-user=freese@uci.edu \
  --mail-type=START,END,FAIL \
  --time=72:00:00" \
  -n
``` -->

Wait for this command to run and make sure the steps that it plans to run are reasonable. After it finishes, run the same command without the `-n` option.

## Pseudochromosome construction

* When you get a new genotype that requires a pseudochromosome (ie has a large humanized insertion), add the relevant metadata to the `pseudochromosome_metadata.tsv` spreadsheet.
* TODO add my flowchart for needs pseudochrom

## Analysis hints

* Visualize your reads from the `data/merge/*.bam` files. BAMs and their indices are stored there. This is very helpful for *verifying incorporation of SNPs from the genotypes*.
* Visualize your coverage tracks (very helpful for interpreting TSSs and TESs) from the `data/merge/*.bw` files.
* TODO find examples on using Swan to visualize transcriptomes in `TODO.ipynb`
