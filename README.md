```bash
snakemake \
  -s Snakefile \
  -j 30 \
  --latency-wait 120 \
  --use-conda \
  -n
  ```

```bash
snakemake -s Snakefile --dag | dot -Tpng > ruledag.png
```

```bash
snakemake \
  -s Snakefile \
  -j 60 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
  ```

  Requirements:

  CLI
    * samtools
    * minimap2
    * faidx
  Python
    * PyRanges
    * TranscriptClean
    * TALON
    * LAPA
    * Cerberus
    * Swan
    * pydeseq2
