#!/bin/bash
#### Working directory. 
WORKDIR="/path/to/your/project"
cd "$WORKDIR"

# /path/to/your/project/trim_damtag is the output directory for trimmed fastq files
# initialize log file
echo "" > trim_damtag/cutadapt_output.txt

# list of FASTQ R1 files (absolute or relative paths)
file_list="/path/to/your/project/raw_fastq_filepaths_damtag.txt"

while read fqfile1; do
    fqfile2=${fqfile1%_R1_001.fastq.gz}_R2_001.fastq.gz
    outfile1="trim_damtag/$(basename "$fqfile1")"
    outfile2="trim_damtag/$(basename "$fqfile2")"


	if [[ ! -e "${outfile1}" ]]
	then 
		# do adapter trimming
		echo "trimming $outfile1"
		cutadapt -a CTGTCTCTTATACACA -A CTGTCTCTTATACACA -m 10 --overlap=1 -o $outfile1 -p $outfile2 $fqfile1 $fqfile2 >> trim_damtag/cutadapt_output.txt

	else 
		echo "skipping existing $outfile1"
	fi
done < "$file_list"
