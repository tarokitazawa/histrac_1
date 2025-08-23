#!/bin/bash
#### Working directory.
WORKDIR="/path/to/your/project"
cd "$WORKDIR"

# /path/to/your/project/trim_dam is the output directory for trimmed fastq files
# initialize log file
echo "" > trim_dam/cutadapt_output.txt

# list of FASTQ R1 files (absolute or relative paths)
file_list="/path/to/your/project/raw_fastq_filepaths_dam.txt"

# Create a temporary directory for intermediate files
mkdir -p temp_trim_dam

while read fqfile1; do
    fqfile2=${fqfile1%_R1_001.fastq.gz}_R2_001.fastq.gz
    outfile1="trim_dam/filtered1_$(basename "$fqfile1")"
    outfile2="trim_dam/filtered1_$(basename "$fqfile2")"
    temp_outfile1="temp_trim_dam/$(basename "$fqfile1")"
    temp_outfile2="temp_trim_dam/$(basename "$fqfile2")"


	if [[ ! -e "${outfile1}.gz" ]]
	then 
		# Filter for 5' TC only and then trim adapters
		echo "filtering for 5' TC and trimming adapters for $outfile1"
		cutadapt -g ^TC -G ^TC -m 10 --discard-untrimmed -o $temp_outfile1 -p $temp_outfile2 $fqfile1 $fqfile2 >> trim_dam/cutadapt_output_temp.txt

		# Check if the temporary files were created
		if [[ -e "$temp_outfile1" && -e "$temp_outfile2" ]]; then
			# Trim adapters from the filtered reads
			cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 10 --overlap=1 -o $outfile1 -p $outfile2 $temp_outfile1 $temp_outfile2 >> trim_dam/cutadapt_output.txt

			# remove temporary files
			rm $temp_outfile1 $temp_outfile2
		else
			echo "Error: Temporary files were not created for $outfile1. Skipping this pair."
		fi
	else 
		echo "skipping existing $outfile1"
	fi
done < "$file_list"
