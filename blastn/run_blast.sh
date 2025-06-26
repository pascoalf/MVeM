# Change path as necessary
blastn -db ./db/mito -query ./ASV.fasta -out ../output/blast_results -outfmt "6 delim=, qacc qlen sseqid sacc slen evalue bitscore score length pident nident mismatch positive gaps staxid ssciname sblastname scomnames"

