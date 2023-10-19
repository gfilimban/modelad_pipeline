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
 snakemake -s Snakefile --forceall --dag | dot -Tpdf > dag.pdf
 snakemake -s Snakefile --forceall --rulegraph | dot -Tpdf > dag.pdf
 ```
