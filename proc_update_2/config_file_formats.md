# Genotype

Genotype metadata is stored in `genotype_metadata.tsv`. If new genotypes are received, they must be put in this file. Below is a description of the necessary columns.
* `genotype`: Name of the genotype. This must match the genotypes listed in `mouse_metadata.tsv` or `analysis_config.tsv`
* `genotype_notes`: Freeform column (not used). Put whatever notes about the genotype here.
* `pseudochrom_needed`: True or False -- whether this genotype has large humanized insertions and thus needs pseudochromosomes added to the reference.
* `pseudochrom`: Comma-separated list of pseudochromosome names. Each unique pseudochromosome listed in this file must have a manually-generated entry of `ref/pseudochrom/<pseudochromosome name>/ref.fa`; such as `ref/pseudochrom/hCLU/ref.fa`.
