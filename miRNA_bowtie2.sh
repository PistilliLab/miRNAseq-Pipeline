#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <fastq_dir> <index> <adapter_sequence> <length_cutoff> <total_threads>"
    exit 1
fi

# Assign input arguments to variables
fastq_dir=$1
index=$2
adapter_sequence=$3
length_cutoff=$4
total_threads=$5

# Create necessary directories
# mkdir -p "${fastq_dir}/fastqc_raw"
# mkdir -p "${fastq_dir}/cutadapt_logs"
# mkdir -p "${fastq_dir}/noadapt_fastqc"
# mkdir -p "${fastq_dir}/noadapt_fastqs"
# mkdir -p "${fastq_dir}/bam_files"
# mkdir -p "${fastq_dir}/bowtie_logs"
mkdir -p "${fastq_dir}/counts_file"
mkdir -p "${fastq_dir}/samtools_logs"
# mkdir -p "${fastq_dir}/log_files"

# # Create a timestamped logfile
# log_timestamp=$(date "+%Y-%m-%d_%H-%M-%S")
# log_command_file="${fastq_dir}/log_files/command_log_${log_timestamp}.log"
# echo "Command executed: $0 $@" > "$log_command_file"

# # Run FastQC on raw sequencing files
# fastqc -o "${fastq_dir}/fastqc_raw" -t $total_threads "${fastq_dir}"/*.fastq

# # Run MultiQC on raw FastQC results
# cd "${fastq_dir}/fastqc_raw"
# multiqc ./
# cd "../"

# echo "Processing fastq files..."

# # Loop through each .fastq file in the directory and run Cutadapt
# num_fastq=$(ls -1 "${fastq_dir}"/*.fastq | wc -l)
# cutadapt_threads=$((total_threads / num_fastq))
# if [ "$cutadapt_threads" -lt 1 ]; then
#     cutadapt_threads=1
# fi

# for input_file in "${fastq_dir}"/*.fastq; do
#     if [ -e "$input_file" ]; then
#         noadapt_file="${input_file%.fastq}_noadapt.fastq"
#         log_file="${fastq_dir}/cutadapt_logs/$(basename "$input_file").log"

#         echo "Processing $input_file with Cutadapt..."
#         cutadapt -j "$cutadapt_threads" -a "$adapter_sequence" -l "$length_cutoff" -o "$noadapt_file" "$input_file" > "$log_file" 2>&1 &

#         # Limit concurrent processes to total_threads
#         while [ $(jobs | wc -l) -ge $total_threads ]; do
#             sleep 1
#         done
#     fi
# done

# echo "Waiting for Cutadapt to complete..."
# wait

# echo "Running FastQC on trimmed files..."
# fastqc -o "${fastq_dir}/noadapt_fastqc" -t $total_threads "${fastq_dir}"/*_noadapt.fastq

# # Run MultiQC on trimmed FastQC results
# cd "${fastq_dir}/noadapt_fastqc"
# multiqc ./
# cd "../"

# # Move noadapt fastqs to their dedicated directory
# mv "${fastq_dir}"/*noadapt.fastq "${fastq_dir}/noadapt_fastqs"

# echo "Performing Bowtie2 alignments..."

alignment_summary_file="${fastq_dir}/bowtie_logs/alignment_summary.tsv"
echo -e "File\tTotal Reads\tReads that Aligned 0 Times\tReads that Aligned Exactly 1 Time\tReads that Aligned >1 Times\tPercent Aligned 0 Times (%)\tPercent Aligned Exactly 1 Time (%)\tPercent Aligned >1 Times (%)\tOverall Alignment Rate (%)" > "$alignment_summary_file"

# Perform alignment to mature miRNAs using Bowtie2
for input_file in "${fastq_dir}/noadapt_fastqs"/*.fastq; do
    if [ -e "$input_file" ]; then
        base_name=$(basename "$input_file" .fastq)
        # sam_file="${fastq_dir}/bam_files/${base_name}.sam"
        # bam_file="${fastq_dir}/bam_files/${base_name}.bam"
        log_file="${fastq_dir}/bowtie_logs/${base_name}.log"
        # echo "Aligning $input_file with Bowtie2..."
        # bowtie2 -p "$total_threads" -q --local --very-sensitive-local -N 0 -L 8 \
        #     --score-min L,4,1 -x "$index" -U "$input_file" -S "$sam_file" > "$log_file" 2>&1

        # echo "Converting SAM to BAM for $sam_file..."
        # samtools view -@ "$total_threads" -bS "$sam_file" | samtools sort -@ "$total_threads" -o "$bam_file"
        # rm "$sam_file"  # Remove the SAM file to save space

        # Extract alignment summary from log file
        total_reads=$(grep "reads; of these:" "$log_file" | awk '{print $1}')
        aligned_0=$(grep "aligned 0 times" "$log_file" | awk '{print $1}')
        aligned_1=$(grep "aligned exactly 1 time" "$log_file" | awk '{print $1}')
        aligned_more=$(grep "aligned >1 times" "$log_file" | awk '{print $1}')
        alignment_rate=$(grep "overall alignment rate" "$log_file" | awk '{print $1}')
        percent_0=$(awk "BEGIN {print ($aligned_0/$total_reads)*100}")
        percent_1=$(awk "BEGIN {print ($aligned_1/$total_reads)*100}")
        percent_more=$(awk "BEGIN {print ($aligned_more/$total_reads)*100}")

        echo -e "$base_name\t$total_reads\t$aligned_0\t$aligned_1\t$aligned_more\t$percent_0\t$percent_1\t$percent_more\t$alignment_rate" >> "$alignment_summary_file"
    fi
done

echo "Generating read count matrix..."
counts_file="${fastq_dir}/counts_file/counts.tsv"

# Initialize the counts file with the header
printf "miRNA" > "$counts_file"
for bam_file in "${fastq_dir}/bam_files"/*.bam; do
    printf "\t%s" "$(basename "$bam_file" .bam)" >> "$counts_file"
done
printf "\n" >> "$counts_file"

# Collect counts for each BAM file in one call
declare -A counts
for bam_file in "${fastq_dir}/bam_files"/*.bam; do
    echo "Processing $bam_file with samtools idxstats..."
    samtools idxstats -@ "$total_threads" "$bam_file" &> "${fastq_dir}/samtools_logs/$(basename "$bam_file" .bam).log"
    while read -r miRNA _ mapped _; do
        counts["$miRNA","$bam_file"]="$mapped"
    done < "${fastq_dir}/samtools_logs/$(basename "$bam_file" .bam).log"
done

# Extract unique miRNA names
first_bam=$(find "${fastq_dir}/bam_files" -maxdepth 1 -type f -name "*.bam" | head -n 1)
miRNAs=$(samtools idxstats -@ "$total_threads" "$first_bam" 2> /dev/null | cut -f1)

# Populate counts table
for miRNA in $miRNAs; do
    printf "%s" "$miRNA" >> "$counts_file"
    for bam_file in "${fastq_dir}/bam_files"/*.bam; do
        printf "\t%s" "${counts["$miRNA","$bam_file"]:-0}" >> "$counts_file"
    done
    printf "\n" >> "$counts_file"
done

echo "Read count matrix saved to $counts_file"

echo "Alignment and read counting completed."
