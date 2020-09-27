# homologr v0.0.1

## New functions:

- adding `diamond_best_hits()`: Perform a DIAMOND best hit search
- adding `diamond_reciprocal_best_hits()`: Perform a DIAMOND best reciprocal hit search
- adding `filter_best_hits()`: Helper function to select best DIAMOND hit based on minimum evalue
- adding `codon_aln()`: Codon Alignments using PAL2NAL
- adding `dnds()`: Genome-wide pairwise inference of synonymous versus non-synonymous substitution rates (dnds)
- adding `dnds_across_multiple_species()`: Genome-wide pairwise inference of synonymous versus non-synonymous substitution rates (dnds) across multiple subject organisms
- adding `import_cds()`: Import the coding sequences from a fasta file
- adding `import_dnds_across_multiple_species()`: Import 'dnds' tables by gene locus and splice varaint for a set of species
- adding `import_proteome()`: Import protein sequences from a 'fasta' file
- adding `internal_dnds()`: Internal dnds computations
- adding `pairwise_aln()`: Pairwise Global Sequence Alignments
- adding `substitutionrate()`: Internal function for dNdS computations
- adding `dnds_across_multiple_species_plot()`: Plot pairwise dN, dS, or dN/dS distributions of orthologous genes across species
