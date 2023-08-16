#!/bin/bash
#SBATCH --job-name=dumb_talon
#SBATCH -p standard
#SBATCH -A seyedam_lab
#SBATCH -c 16
#SBATCH --error=slurm-%J.err
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mem 360GB
#SBATCH --time=24:00:00

sam=/dfs6/pub/freese/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/5xBIN1_HO_F_4_months_HC_3_labeled_merged.bam

touch test_config.csv
rm test_config.csv
touch test_config.csv
echo "test,test,pacbio,"${sam} >> test_config.csv

ref_gtf=/dfs6/pub/freese/mortazavi_lab/bin/modelad_pipeline/data/230516/cerberus/ca_vM21.gtf
talon_initialize_database \
  --f ${ref_gtf} \
  --g mm10 \
  --a ca_vM21 \
  --o test

talon \
  --f test_config.csv \
  --db test.db \
  --build mm10 \
  -t 16 \
  --o test_2
