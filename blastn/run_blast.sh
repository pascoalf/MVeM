# Change path as necessary
blastn -db nt -query ASV_eDNA.fasta -out blast_results -outfmt "6 delim=, qacc qlen sseqid sacc slen evalue bitscore score length pident nident mismatch positive gaps staxid ssciname sblastname scomnames skingdoms" -remote

