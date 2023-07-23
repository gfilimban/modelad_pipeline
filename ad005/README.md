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
