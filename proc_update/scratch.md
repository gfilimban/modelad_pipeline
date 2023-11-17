```bash
fastq=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/human_ref/pseudochrom/hCLU/CLU/CLU_t.fastq
fasta=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/human_ref/pseudochrom/hCLU/CLU/CLU_t.fasta
# https://www.biostars.org/p/85929/
sed -n '1~4s/^@/>/p;2~4p' $fastq > $fasta

ref=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/ref/pseudochrom/hCLU/ref.fa
out=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/human_ref/pseudochrom/hCLU/CLU/CLU_blat.psl
blat $ref $fasta -t=dna -q=dna -out=psl $out

bed=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/human_ref/pseudochrom/hCLU/CLU/CLU_blat.bed
pslToBed $out $bed

gp=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/human_ref/pseudochrom/hCLU/CLU/CLU_blat.gp
bedToGenePred $bed $gp

gtf=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/human_ref/pseudochrom/hCLU/CLU/CLU_blat.gtf
genePredToGtf $gp $gtf

# conda activate igvtools
# sort_gtf=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/human_ref/pseudochrom/hCLU/CLU/CLU_blat_sorted.gtf
# igvtools $gtf $sort_gtf
```

```python
import glob
import pyranges as pr
import shutil
r = 'data/talon/*/*gtf'
# for f in glob.glob(r):
#   shutil.copyfile(f, f+'_back')

for f in glob.glob(r):
  df = pr.read_gtf(f, rename_attr=True).df
  df.drop('source_attr', axis=1, inplace=True)
  pr.PyRanges(df).to_gtf(f)


```

```bash
lapa_correct_talon \
  --links data/lapa/ad003/5xCLU-h2kbKI_HO_F_4_months_HC/tss_to_tes_links.csv \
  --read_annot data/talon/ad003/5xCLU-h2kbKI_HO_F_4_months_HC_talon_read_annot.tsv \
  --gtf_input data/talon/ad003/5xCLU-h2kbKI_HO_F_4_months_HC_talon.gtf \
  --gtf_output data/lapa/ad003/5xCLU-h2kbKI_HO_F_4_months_HC/talon_corrected/corrected.gtf \
  --abundance_input data/talon/ad003/5xCLU-h2kbKI_HO_F_4_months_HC_talon_abundance_filtered.tsv \
  --abundance_output data/lapa/ad003/5xCLU-h2kbKI_HO_F_4_months_HC/talon_corrected/filtered_abundance_corrected.tsv \
  --keep_unsupported
 ```
