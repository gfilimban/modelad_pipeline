If preliminary / exploratory analysis was needed for any study / experiment, then there will be materials in the corresponding folder. Ie
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
--cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n

snakemake \
-s 231022_ad006_snakefile.smk \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
```
