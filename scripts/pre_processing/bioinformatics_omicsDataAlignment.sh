#!/bin/bash

# BASH shell script to align omics data for biostatistical analysis

# move to the directory with the data
#cd "/PATH/TO/DATA"

# convert gff3 file to gtf
gffread -E -F -T Tribolium_castaneum.gff3 -o Tribolium.gtf

# generate reference genome build files
hisat2-build Tribolium_castaneum.genome.fa TriboliumBuild

# align trimmed samples to the refence genome
# cntrl samples 4h
hisat2 -q -x TriboliumBuild -U SRR8288561.trimmed.fq.gz -S SRR8288561_accepted_hits.sam --summary-file SRR8288561_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288562.trimmed.fq.gz -S SRR8288562_accepted_hits.sam --summary-file SRR8288562_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288563.trimmed.fq.gz -S SRR8288563_accepted_hits.sam --summary-file SRR8288563_alignedSummary.txt
# cntrl samples 24h
hisat2 -q -x TriboliumBuild -U SRR8288558.trimmed.fq.gz -S SRR8288558_accepted_hits.sam --summary-file SRR8288558_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288567.trimmed.fq.gz -S SRR8288567_accepted_hits.sam --summary-file SRR8288567_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288568.trimmed.fq.gz -S SRR8288568_accepted_hits.sam --summary-file SRR8288568_alignedSummary.txt
# treat samples 4h
hisat2 -q -x TriboliumBuild -U SRR8288564.trimmed.fq.gz -S SRR8288564_accepted_hits.sam --summary-file SRR8288564_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288557.trimmed.fq.gz -S SRR8288557_accepted_hits.sam --summary-file SRR8288557_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288560.trimmed.fq.gz -S SRR8288560_accepted_hits.sam --summary-file SRR8288560_alignedSummary.txt
# treat samples 24h
hisat2 -q -x TriboliumBuild -U SRR8288559.trimmed.fq.gz -S SRR8288559_accepted_hits.sam --summary-file SRR8288559_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288565.trimmed.fq.gz -S SRR8288565_accepted_hits.sam --summary-file SRR8288565_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288566.trimmed.fq.gz -S SRR8288566_accepted_hits.sam --summary-file SRR8288566_alignedSummary.txt
