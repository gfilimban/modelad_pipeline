* Need to add 5xWT and WT from AD003 to SwanGraph

```bash
python make_sg.py
```

* Make sure that the raw reads have the SNP

```bash
cd ~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/
gtf=~/mortazavi_lab/bin/modelad_pipeline/ref/annot.gtf
grep "Bin1" $gtf | head -1
# chr18:32377230-32435737

module load samtools

# sort
bam=5xBIN1_HO_F_4_months_HC_1_labeled_merged.bam
sort_bam=5xBIN1_HO_F_4_months_HC_1_labeled_merged_sorted.bam
mini_sam=5xBIN1_HO_F_4_months_HC_1_labeled_merged_sorted_mini.sam
mini_sam_sorted=5xBIN1_HO_F_4_months_HC_1_labeled_merged_sorted_mini_sorted.sam
mini_bam=5xBIN1_HO_F_4_months_HC_1_labeled_merged_sorted_mini.bam

samtools sort $bam > $sort_bam

# index
samtools index $sort_bam

# get subset, then index again
samtools view -h $sort_bam chr18 > $mini_sam
samtools sort $mini_sam > $mini_sam_sorted
samtools view -hb $mini_sam_sorted > $mini_bam
samtools index $mini_bam



# sort
bam=BIN1_HO_F_4_months_HC_1_labeled_merged.bam
sort_bam=BIN1_HO_F_4_months_HC_1_labeled_merged_sorted.bam
mini_sam=BIN1_HO_F_4_months_HC_1_labeled_merged_sorted_mini.sam
mini_sam_sorted=BIN1_HO_F_4_months_HC_1_labeled_merged_sorted_mini_sorted.sam
mini_bam=BIN1_HO_F_4_months_HC_1_labeled_merged_sorted_mini.bam

samtools sort $bam > $sort_bam

# index
samtools index $sort_bam

# get subset, then index again
samtools view -h $sort_bam chr18 > $mini_sam
samtools sort $mini_sam > $mini_sam_sorted
samtools view -hb $mini_sam_sorted > $mini_bam
samtools index $mini_bam

```
<!-- chr18:32377230-32435737 -->








## 230815 compare 3' end usage in short and long read

```bash
# get bigwigs for each sr bam file
d=/share/crsp/lab/model-ad/share/bulkRNA/5x_GWAS/5xBin1/5xBin1_4mo/


```

```bash
snakemake \
  -s Snakefile \
  -j 30 \
  --latency-wait 120 \
  --use-conda \
  -n
  ```

```bash
conda activate viz_snakemake
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


FUSION READS, WHAT ARE YOU DOING HERE

```bash
module load samtools
chr11:102881298-102886826

cd ~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon
for f in *BIN1*merged.bam
do
  # samtools sort --threads 8 -O bam $f > ${f}_sorted.bam
  # samtools index ${f}_sorted.bam
  samtools view -hb ${f}_sorted.bam chr11:102881298-102886826 > ${f}_subset.bam
  samtools sort --threads 8 -O bam ${f}_subset.bam > ${f}_subset_sorted.bam
  samtools index ${f}_subset_sorted.bam
done
```

# get just that problem read
```bash
module load samtools

bam=~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/5xBIN1_HO_F_4_months_HC_3_labeled_merged.bam_subset_sorted.bam

read_name=212bf012-20aa-4005-a74b-b3db02810488_1
output_sam=~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/test_subset.sam
samtools view -HS $bam > $output_sam
samtools view -S $bam | grep $read_name >> $output_sam

```

# run a stupid little talon run
```bash
sam=~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/test_subset.sam

touch test_config.csv
echo "test,test,pacbio,"${sam} >> test_config.csv

ref_gtf=/Users/fairliereese/Documents/programming/mortazavi_lab/bin/modelad_pipeline/data/230516/cerberus/ca_vM21.gtf
talon_initialize_database \
  --f ${ref_gtf} \
  --g mm10 \
  --a ca_vM21 \
  --o test

talon \
  --f test_config.csv \
  --db test.db \
  --build mm10 \
  -t 1 \
  --o test_2
```

# run another stupid little talon run on the whole file with the new settings
```bash
sam=~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/5xBIN1_HO_F_4_months_HC_3_labeled_merged.bam

touch test_config.csv
echo "test,test,pacbio,"${sam} >> test_config.csv

ref_gtf=/Users/fairliereese/Documents/programming/mortazavi_lab/bin/modelad_pipeline/data/230516/cerberus/ca_vM21.gtf
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

talon_filter_transcripts \
    --db test.db \
    -a ca_vM21 \
    --maxFracA=0.5 \
    --minCount=5 \
    --o output_list.csv

talon_create_GTF \
    --db test.db \
    -b mm10 \
    -a ca_vM21 \
    --whitelist output_list.csv \
    --observed \
    --o test_2_gtf

```

```bash
# extract only the reads that were hits to the dumb new gene that talon created
module load samtools
bam=~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/5xBIN1_HO_F_4_months_HC_3_labeled_merged.bam
read_ids=gene_id_110958_reads.txt
samtools view -N $read_ids -bh $bam > weird_gene.bam
samtools sort -O bam weird_gene.bam > weird_gene_sorted.bam
samtools index weird_gene_sorted.bam
```

All of the reads that were assigned to the weird gene were monoexonic, with the exception of one read.
Maybe we should consider gene status (ie known or novel) before considering overlap as a a heuristic for calling genes from weird reads, cause this will also affect the gene-level quantification
