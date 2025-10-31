#!/bin/bash
for r1 in *_R1_*.fastq.gz; do
    r2=${r1/_R1_/_R2_}
    cutadapt \
	-a "AGACGAGAAGACCCTATG;e=0.15;o=5...GGGATAACAGCGCAATCC;e=0.15;o=5" \
	-A "GGATTGCGCTGTTATCCC;e=0.15;o=5...CATAGGGTCTTCTCGTCT;e=0.15;o=5" \
	-q 20 \
	-o Trimmed/${r1} \
	-p Trimmed/${r2} ${r1} ${r2}
done
