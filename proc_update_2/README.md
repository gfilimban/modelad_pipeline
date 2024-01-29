# MODEL-AD Long-read RNA-seq pipeline

Preparing your processing environment:
* TranscriptClean
* gencode_utr_fix
* igvtools
* deeptools

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

<!-- If preliminary / exploratory analysis was needed for any study / experiment, then there will be materials in the corresponding folder. Ie
* ad008 folder, to check out quickly from mapped reads worked

```bash
snakemake \
  -s Snakefile \
  -j 30 \
  --latency-wait 120 \
  --use-conda \
  -n
  ```

  ```bash
 conda activate snakemake_vis
 snakemake -s Snakefile --forceall --dag | dot -Tpdf > pipeline_dag.pdf
 snakemake -s Snakefile --forceall --rulegraph | dot -Tpdf > pipeline_rulegraph.pdf
 ```

```bash
snakemake \
-s Snakefile \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A model-ad_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n

snakemake \
-s 231122_hapoe4_snakefile.smk \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n

snakemake \
-s 231022_ad006_snakefile.smk \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
```

```bash
talon_filter_transcripts \
        --db hTREM2-R47HKI_HET_F_8_weeks_HC_talon.db \
        -a vM21 \
        --maxFracA 0.5\
        --minCount 5 \
        --minDatasets 2 \
        --filter_known \
        --o temp
``` -->
