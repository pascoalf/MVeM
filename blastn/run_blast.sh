# Change path as necessary
blastn -db nt -query ASV.fasta -out blast_results -outfmt "6 qacc qlen sseqid sacc slen evalue bitscore score length pident nident mismatch positive gaps staxid ssciname sblastname scomnames skingdoms" -evalue 1e-05 -perc_identity 99 -qcov_hsp_perc 80 -remote
