#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <fastq_dir> <sample_list> <index_dir> <hairpin_fasta> <mature_fasta> <species> <total_threads>"
    exit 1
fi

# Assign input arguments to variables
fastq_dir=$1
sample_list=$2
index_dir=$3
hairpin_fasta=$4
mature_fasta=$5
species=$6
total_threads=$7

source ~/.bashrc

# Create a directory for raw_fastqc
mkdir "${fastq_dir}/fastqc_raw"

# Run fastqc on raw sequencing files
fastqc -o "${fastq_dir}/fastqc_raw" -t $total_threads "${fastq_dir}"/*.fastq

# Run multiqc on all raw fastqc files
cd "${fastq_dir}/fastqc_raw"
multiqc ./
cd "../"

# Create a directory for log files
mkdir "${fastq_dir}/cutadapt_logs"

echo "Processing fastq files..."

# Loop through each .fastq file in the directory
for input_file in "$fastq_dir"/*.fastq; do
    if [ -e "$input_file" ]; then
        noadapt_file="${input_file%.fastq}_noadapt.fastq"
        log_file="${fastq_dir}/cutadapt_logs/$(basename "$input_file").log" # Create a log file path

        echo "Processing $input_file..."

        cutadapt -j 2 -a TGGAATTCTCGGGTGCCAAGG -m 26 -M 36 -o "$noadapt_file" "$input_file" > "$log_file" 2>&1 &

        # Keep track of running processes
        while [ $(jobs | wc -l) -ge $total_threads ]; do
            sleep 1
        done
    fi
done

echo "Waiting for cutadapt to complete..."

# Wait for all processes to finish
wait

# Create a directory for log files
mkdir "${fastq_dir}/trimmed_logs"

echo "Trimming fastq files..."

# Loop through each .fastq file in the directory
for input_file in "$fastq_dir"/*_noadapt.fastq; do
    if [ -e "$input_file" ]; then
        trimmed_file="${input_file%.fastq}_trimmed.fastq"
        log_file="${fastq_dir}/trimmed_logs/$(basename "$input_file").log" # Create a log file path

        echo "Processing $input_file..."

        cutadapt -j 2 -u 5 -u -6 -o "$trimmed_file" "$input_file" > "$log_file" 2>&1 &

        # Keep track of running processes
        while [ $(jobs | wc -l) -ge $total_threads ]; do
            sleep 1
        done
    fi
done

echo "Waiting for trimming to complete..."

# Wait for all processes to finish
wait

# Run fastqc on all trimmed fastqs
mkdir "${fastq_dir}/trimmed_fastqc"
fastqc -o "${fastq_dir}/trimmed_fastqc" -t $total_threads "${fastq_dir}"/*_trimmed.fastq

echo "FastQC completed."

# Run multiqc on all trimmed fastqc files
cd "${fastq_dir}/trimmed_fastqc"
multiqc "./"
cd "../"

echo "MultiQC completed."

# Move noadapt fastqs to a new directory
mkdir "${fastq_dir}/noadapt_fastqs"
mv "${fastq_dir}"/*noadapt.fastq "${fastq_dir}/noadapt_fastqs"

# Move trimmed fastqs to a new directory
mkdir "${fastq_dir}/trimmed_fastqs"
mv "${fastq_dir}"/*_trimmed.fastq "${fastq_dir}/trimmed_fastqs"

# Running mapper.pl
echo "Running mapper.pl..."

# Create a directory for collapsed files, arf files, and mapper.pl logs
mkdir "${fastq_dir}/collapsed_fasta"
mkdir "${fastq_dir}/arf_file"
mkdir "${fastq_dir}/mapper_logs"

# Extract the genome name from the index directory
genome_name=$(basename "$index_dir"/*.1.ebwt)
genome_name=${index_dir}/${genome_name%.1.ebwt}

# Move to index directory
cd "${fastq_dir}/trimmed_fastqs"

# Run mapper.pl with sample list
mapper.pl "$sample_list" -d -e -h -r 10000 -m -j -l 14 -p "$genome_name" -s mapper_collapsed.fa -t mapper_aligned.arf -o "$total_threads" -v > "${fastq_dir}/mapper_logs/mapper.log" 2>&1

# Move files to respective directories
mv mapper_collapsed.fa "${fastq_dir}/collapsed_fasta"
mv mapper_aligned.arf "${fastq_dir}/arf_file"
mv mapper_logs "${fastq_dir}/mapper_logs"

# Delete bowtie.log file in index directory
rm bowtie.log

# Return to the fastq directory
cd "${fastq_dir}"

echo "mapper.pl completed."

# Make a quantifier and quantifier_logs directories
mkdir "${fastq_dir}/quantifier"
mkdir "${fastq_dir}/quantifier/quantifier_logs"
cd "${fastq_dir}/quantifier"

echo "Running quantifier.pl..."

# Run quantifier.pl
quantifier.pl -d -p ${hairpin_fasta} -m ${mature_fasta} -r ${fastq_dir}/collapsed_fasta/mapper_collapsed.fa -t ${species} > "${fastq_dir}/quantifier/quantifier_logs/quantifier.log" 2>&1

echo "quantifier.pl completed."

echo "Processing completed."
