# MODEL-AD Long-read RNA-seq pipeline

## Getting started

* First, download this repository
```bash
git clone git@github.com:fairliereese/modelad_pipeline.git
```
* Build your necessary environments / download necessary software: [Link](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/environment.md).

## Preprocessing

Best practices:
* Get familiar with using `tmux` to stash your sessions while performing preprocessing
* Always use `srun` on HPC3 before running Snakemake or you'll get yelled at
* Make sure to run `conda activate snakemake` before running any Snakemake commands

Checklist before running:
* When you get new reads, add the paths to FASTQ files to [`config.tsv`](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/config.tsv).
* Add new mice / tissues to [`mouse_metadata.tsv`](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/mouse_metadata.tsv)
* Add new genotypes to [`genotype_metadata.tsv`](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/genotype_metadata.tsv)
* Add new pseudochromosomes (if needed) to [`pseudochromosome_metadata.tsv`](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/pseudochromosome_metadata.tsv)
* Make sure to save, commit, and push all your changes to metadata files back to this repo.
* Navigate to the public pipeline directory on MODEL-AD CRSP (`/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update_2/`) and run `git pull` to add metadata file updates.


Run the preprocessing pipeline with the following command (replace `{your_email}` with your email you want updates at!). `Snakefile` will run the preprocessing pipeline through LAPA (mapping->TranscriptClean->TALON->LAPA):
```bash
conda activate snakemake
snakemake \
-s Snakefile \
-j 200 \
--latency-wait 120 \
--use-conda \
--rerun-triggers mtime \
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
-s Snakefile \
-j 200 \
--latency-wait 120 \
--use-conda \
  -n

snakemake \
-s Snakefile \
-j 200 \
--latency-wait 120 \
--use-conda \
--rerun-triggers mtime \
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

#### Analyzing preprocessing output files

The analysis goal of the preprocessing output is to confirm that:
* *Splicing of mutated / altered genes has not been changed*
* *SNPs of altered genes have been incorporated*.
* Visualize your reads from the `/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update_2/data/merge/*.bam` files. BAMs and their indices are stored there.
* For genotypes with pseudochromosomes, also download the genome reference (`*.fa` and `*.fai`) and transcriptome reference (`*.gtf`) for your target genotype (in this folder: `/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update_2/ref/genome/<your genotype>/`)
* I recommend using [IGV](https://www.igv.org/) to visualize your reads.

* Visualize your coverage tracks (very helpful for interpreting TSSs and TESs) from the `data/merge/*.bw` files.

## Downstream analysis

For the majority of you, you will only be running this part of the pipeline. This part of the pipeline will run some basic analyses on the long-read data for you including differential gene and transcript expression tests and differential isoform usage tests.

* First, make sure you _fork_ and _clone (`git clone`) this repo_ so that you have access to all the relevant metadata and you can edit the `analysis_config.tsv`.
* When you want to run a new analysis, edit or create a new [`analysis_config.tsv`](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/analysis_config.tsv). Check the [`mouse_metadata.tsv`](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/mouse_metadata.tsv) file for the relevant genotypes and studies to include in your table.
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
--ignore-incomplete \
  -n

conda activate modelad_snakemake_2
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
   --ignore-incomplete -n

   conda activate modelad_snakemake_2
   snakemake \
   -s 240211_nm_snakefile_analysis.smk \
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

## Visualize the pipelines execution

```bash
conda activate snakemake_vis

# preprocessing
snakemake -s Snakefile --forceall --dag | dot -Tpdf > pipeline_dag.pdf # entire DAG
snakemake -s Snakefile --forceall --rulegraph | dot -Tpdf > pipeline_rulegraph.pdf # abbreviated DAG (only 1 box / rule)

# preprocessing
snakemake -s Snakefile --forceall --dag | dot -Tpdf > analysis_pipeline_dag.pdf # entire DAG
snakemake -s Snakefile --forceall --rulegraph | dot -Tpdf > analysis_pipeline_rulegraph.pdf # abbreviated DAG (only 1 box / rule)
```

## Pseudochromosome construction

* When you get a new genotype that requires a pseudochromosome (ie has a large humanized insertion), add the relevant metadata to the [`pseudochromosome_metadata.tsv`](https://github.com/fairliereese/modelad_pipeline/blob/main/proc_update_2/pseudochromosome_metadata.tsv) spreadsheet.
* TODO add my flowchart for needs pseudochrom

## Analysis hints

* TODO find examples on using Swan to visualize transcriptomes in `TODO.ipynb`
