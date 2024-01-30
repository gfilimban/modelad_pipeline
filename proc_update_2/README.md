# MODEL-AD Long-read RNA-seq pipeline


## Preprocessing
Preparing your processing environment:
* TranscriptClean
* gencode_utr_fix
* igvtools
* deeptools
* TODO add yaml config files
* TODO mention qi and tmux


* When you get new reads, add the paths to FASTQ files to `config.tsv`.
* `Snakefile` will run the preprocessing pipeline through LAPA.

Run the preprocessing pipeline with the following command:
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
