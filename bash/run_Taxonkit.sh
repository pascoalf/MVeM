cat blast_results | taxonkit reformat2 -I 15 -r "Unassigned" -f "{domain|acellular root|superkingdom}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}" | tee blast_tax_results
