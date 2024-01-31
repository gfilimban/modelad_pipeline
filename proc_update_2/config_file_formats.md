# Preprocessing


## Genotype

Genotype metadata is stored in `genotype_metadata.tsv`. If new genotypes are received, they must be put in this file. Below is a description of the necessary columns.
* `genotype`: Name of the genotype. This must match the genotypes listed in `mouse_metadata.tsv` or `analysis_config.tsv`
* `genotype_notes`: Freeform column (not used). Put whatever notes about the genotype here.
* `pseudochrom_needed`: True or False -- whether this genotype has large humanized insertions and thus needs pseudochromosomes added to the reference.
* `pseudochrom`: Comma-separated list of pseudochromosome names. Each unique pseudochromosome listed in this file must have a manually-generated entry of `ref/pseudochrom/<pseudochromosome name>/ref.fa`; such as `ref/pseudochrom/hCLU/ref.fa`.

# Analysis

## Analysis config

When you start a new project and want to analyze a group of mice / genotypes together, edit your own `analysis_config.tsv` file. Below is a description of the necessary columns:
* `analysis`: Name of the analysis you're running. Will dictate where the output files are stored.
*  `genotype`: Name of genotype to include. This must match genotypes listed in `genotype_metadata.tsv` and `mouse_metadata.tsv`.
* `study`: Name of the study to pull the corresponding mice of `genotype` from. This is useful if there is more than one set of mice of the same genotype, but you only want to include the one. Refer to the studies listed in the `Experiment` column of [this spreadsheet](https://docs.google.com/spreadsheets/d/1jR0lWpx3t42vlJfhsdu5OAZRJ4DnRhWUdAJNXhMZzmU/edit?usp=sharing).
* `genotype_alias`: (Optional) A new name to be used to display and group your datasets for comparisons. For instance, if you have two studies included that both have WTs but the `genotype`s don't match (ie "hABKI-Swe-WT" and "hTREM2KI-WT", which are both just Black 6, but do not match based on their original genotype names). If these are not filled in the original genotype from the `genotype` column will be used.
