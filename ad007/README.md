* Need to add 5xWT and WT from AD003 to SwanGraph

```bash
python make_sg.py
```

* Make sure that the raw reads have the SNP

```bash
cd ~/mortazavi_lab/bin/modelad_pipeline/data/230516/talon/
gtf=~/mortazavi_lab/bin/modelad_pipeline/ref/annot.gtf
grep "Abca7" $gtf | head -1
# chr10:79996494-80015572

module load samtools

# sort
bam=ABCA7_HO_F_4_months_HC_2_labeled_merged.bam
sort_bam=ABCA7_HO_F_4_months_HC_2_labeled_merged_sorted.bam
mini_sam=ABCA7_HO_F_4_months_HC_2_labeled_merged_sorted_mini.sam
mini_sam_sorted=ABCA7_HO_F_4_months_HC_2_labeled_merged_sorted_mini_sorted.sam
mini_bam=ABCA7_HO_F_4_months_HC_2_labeled_merged_sorted_mini.bam

samtools sort $bam > $sort_bam

# index
samtools index $sort_bam

# get subset, then index again
samtools view -h $sort_bam chr10 > $mini_sam
samtools sort $mini_sam > $mini_sam_sorted
samtools view -hb $mini_sam_sorted > $mini_bam
samtools index $mini_bam


```
