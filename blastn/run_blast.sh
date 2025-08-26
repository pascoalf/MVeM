# Change path as necessary
blastn -db nt -query ASV_atlantida.fasta -out blast_results_atlantida -outfmt "6 delim=, qacc qlen sseqid sacc slen evalue bitscore score length pident nident mismatch positive gaps staxid ssciname sblastname scomnames skingdoms" -evalue 1e-05 -perc_identity 99 -qcov_hsp_perc 80 -remote
