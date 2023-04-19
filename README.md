```bash
snakemake \
  -s Snakefile \
  -j 10 \
  --latency-wait 120 \
  -n
  ```


```bash
snakemake \
  -s snakemake/human_snakefile.smk \
  -j 10 \
  --latency-wait 120 \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END, --time=72:00:00" -n
  ```
