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
