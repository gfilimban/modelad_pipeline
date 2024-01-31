# Preprocessing

## Genotype metadata

Genotype metadata is stored in `genotype_metadata.tsv`. If new genotypes are received, they must be put in this file. Below is a description of the necessary columns.
* `genotype`: Name of the genotype. This must match the genotypes listed in `mouse_metadata.tsv` or `analysis_config.tsv`
* `genotype_notes`: Freeform column (not used). Put whatever notes about the genotype here.
* `pseudochrom_needed`: True or False -- whether this genotype has large humanized insertions and thus needs pseudochromosomes added to the reference.
* `pseudochrom`: Comma-separated list of pseudochromosome names. Each unique pseudochromosome listed in this file must have a manually-generated entry of `ref/pseudochrom/<pseudochromosome name>/ref.fa`, such as `ref/pseudochrom/hCLU/ref.fa`, and must match entries in `pseudochromosome_metadata.tsv`.

## Pseudochromosome metadata

Metadata about each constructed pseudochromosome is stored in `pseudochromosome_metadata.tsv`. When new pseudochromosomes are constructed, their sequences must be stored in `ref/pseudochrom/<pseudochromosome name>/ref.fa`, such as `ref/pseudochrom/hCLU/ref.fa`, and must be added to this table, as well as [the table on OneDrive](https://ucirvine.sharepoint.com/:x:/r/sites/biosci-ad-model/Shared%20Documents/BDMC/HUMANIZED%20ALLELES%20INFO/humanized_allele_refs.xlsx?d=wa587942b349e4826840aaa58463bc713&csf=1&web=1&e=ZF8JYd). Information on what the content of the pseudochromosomes consist of are in [this OneDrive folder](https://ucirvine.sharepoint.com/:f:/s/biosci-ad-model/Eo2rONqobYhJt8BzpHpJ1xEBHV-q8HOoCYmX7UPeeZI-zQ?e=g6HweL). Below is a description of the necessary columns.
* `pseudochrom`: Name of the pseudochromosome. Must match those listed in `genotype_metadata.tsv` as well as the chromosome name given in `ref/pseudochrom/<pseudochromosome name>/ref.fa`.
* `human_gene`: Name of human gene, if any, that is present in this pseudochromosome.
* `mouse_gene`: Name of mouse gene, if any, that is present in this pseudochromosome.
* `locus_type`: Type of insertion made. Either none (for the dummy entry); "chimeric", for genes that have some human and some mouse exons; or "human", for genes that are entirely human.
* `notes`: Free text notes regarding the pseudochromosome.

## Mouse metadata

Metadata about each mouse is stored in `mouse_metadata.tsv`. If new mice are received, they must be put in this file. Refer to [this spreadsheet](https://docs.google.com/spreadsheets/d/1jR0lWpx3t42vlJfhsdu5OAZRJ4DnRhWUdAJNXhMZzmU/edit?usp=sharing) for mouse IDs, ages, sexes, etc. Below is a description of the necessary columns.
* `mouse_id`: ID of the mouse.
* `sex`: Sex of the mouse, either "M" or "F".
* `age`: Age of the mouse. Enter these in a regular manner! For example, "4mo" and "4 months" will not be interpreted the same way. Refer to earlier entries in the file for how to format these.
* `tissue`: Tissue from the mouse. Enter these in a regular manner! For example, "Hippocampus" and "HC" will not be interpreted the same way. Refer to earlier entries in the file for how to format these.
* `genotype`: Genotype of the mouse. This must match the genotypes listed in `genotype_metadata.tsv` or `analysis_config.tsv`

## FASTQ config

When you get new data (FASTQ) to process, add entries to the `config.tsv` file. Below is a description of the necessary columns:
* `platform`: The name of the platform that the data was generated with. Currently, just "ONT".
* `fname`: Name of the `.fastq` or `.fastq.gz` file to be processed. NOTE: this pipeline is expecting to pull out specific pieced of information from this file name, including the study ("ad###", ie "ad003"), and mouse id (ie "10721").  These relationships will be insured in a future update when basecalling and trimming are added to this pipeline.

# Analysis

## Analysis config

When you start a new project and want to analyze a group of mice / genotypes together, edit your own `analysis_config.tsv` file. Below is a description of the necessary columns:
* `analysis`: Name of the analysis you're running. Will dictate where the output files are stored.
*  `genotype`: Name of genotype to include. This must match genotypes listed in `genotype_metadata.tsv` and `mouse_metadata.tsv`.
* `study`: Name of the study to pull the corresponding mice of `genotype` from. This is useful if there is more than one set of mice of the same genotype, but you only want to include the one. Refer to the studies listed in the `Experiment` column of [this spreadsheet](https://docs.google.com/spreadsheets/d/1jR0lWpx3t42vlJfhsdu5OAZRJ4DnRhWUdAJNXhMZzmU/edit?usp=sharing).
* `genotype_alias`: (Optional) A new name to be used to display and group your datasets for comparisons. For instance, if you have two studies included that both have WTs but the `genotype`s don't match (ie "hABKI-Swe-WT" and "hTREM2KI-WT", which are both just Black 6, but do not match based on their original genotype names). If these are not filled in the original genotype from the `genotype` column will be used.
