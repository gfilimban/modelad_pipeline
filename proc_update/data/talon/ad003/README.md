```python
import pandas as pd
annots = ['5xCLU-h2kbKI_HO_F_4_months_HC_talon_read_annot.tsv',
          'CLU-h2kbKI_HO_F_4_months_HC_talon_read_annot.tsv']
genotypes = ['5xCLU', 'CLU']
for f, g in zip(annots, genotypes):
    df = pd.read_csv(f, sep='\t')
    temp = df.loc[df.chrom=='hCLU']
    temp[['read_name', 'transcript_novelty']].groupby('transcript_novelty').count()
    temp2 = temp[['gene_ID', 'transcript_ID']].drop_duplicates()
    temp2.to_csv(f'hCLU_{g}_pass_list.csv', index=False, header=False)
```

read_name
transcript_novelty
Antisense                   1
ISM                       892
Known                     177
NIC                      2435
NNC                       175

read_name
transcript_novelty
Antisense                   4
Genomic                     2
ISM                      1217
Known                     245
NIC                      3521
NNC                       220

```bash
conda activate snakemake
talon_create_GTF \
        --db 5xCLU-h2kbKI_HO_F_4_months_HC_talon.db \
        -a vM21 \
        -b mm10 \
        --whitelist 5xCLU_pass_list.csv \
        --o 5xCLU

conda activate igvtools
igvtools sort 5xCLU_talon.gtf 5xCLU_sorted.gtf
igvtools index 5xCLU_sorted.gtf

conda activate snakemake
talon_create_GTF \
        --db CLU-h2kbKI_HO_F_4_months_HC_talon.db \
        -a vM21 \
        -b mm10 \
        --whitelist CLU_pass_list.csv \
        --o CLU

conda activate igvtools
igvtools sort CLU_talon.gtf CLU_sorted.gtf
igvtools index CLU_sorted.gtf
```
