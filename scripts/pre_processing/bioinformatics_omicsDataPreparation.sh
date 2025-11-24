#!/bin/bash

# BASH shell script to prepare omics data for biostatistical analysis

# move to the directory with the data
#cd "/PATH/TO/DATA"

# qc the raw samples
for rawFile in *".fastq.gz"; do
	fastqc $rawFile --extract
done

# aggregate the reports in the current directory
multiqc ./*fastqc* -n "multiqc_raw"

# setup the adapters file path
adapterPath="/afs/crc.nd.edu/x86_64_linux/bio/Trimmomatic/0.32/adapters/TruSeq3-PE.fa"

# setup the phred quality score: 33 for newer Illumina 1.9 data or 64 for older Illumina 1.5 data
# the version of Illumina can be found in the fastqc_data.txt file from the fastqc results directory
score=33

# trim the samples
for sampleFile in *".trimmed.fq.gz"; do
	# remove the file extension
	fileName=$(echo $sampleFile | sed "s/\.fastq\.gz//g")
	# trim a sample based on QC reports 
	# use -threads argument to improve run time, if multi-threading is available on your machine
	trimmomatic SE -phred"$score" $sampleFile $fileName.trimmed.fq.gz ILLUMINACLIP:"$adapterPath":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10
done

# qc the trimmed samples
for trimmedFile in *".trimmed.fq.gz"; do
	fastqc $trimmedFile --extract
done

# aggregate the reports in the current directory
multiqc ./*"trimmed"*"fastqc"* -n "multiqc_trimmed"
