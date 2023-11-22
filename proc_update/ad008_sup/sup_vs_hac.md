## Number of sup trimmed reads
```bash
d=/share/crsp/lab/model-ad/share/LRRNA_nanopore/ad008_hABKISWEIB_8wk/trim_fastq/
for f in ${d}*number_trimmed_reads*txt
do 
    echo $f
    cat $f
    echo ""
done
```

## Number of hac trimmed reads
```bash
d=/share/crsp/lab/model-ad/share/LRRNA_nanopore/ad008_hABKISWEIB_8wk/trim_fastq_hac_new/
for f in ${d}*number_trimmed_reads*txt
do 
    echo $f
    cat $f
    echo ""
done
```

## Number of sup mapped / tc'd reads
```bash
d=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/ad008_sup/ad008_sup/data/merge/ad008/
module load samtools
for f in ${d}*bam
do
    echo $f
    samtools view -c $f
    echo ""
done
```

## Number of sup mapped / tc'd reads
```bash
d=/share/crsp/lab/model-ad/share/freese/modelad_pipeline/proc_update/ad008_hac/ad008_hac/data/merge/ad008/
module load samtools
for f in ${d}*bam
do
    echo $f
    samtools view -c $f
    echo ""
done
```