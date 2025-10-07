#!/bin/bash


###			QUALITY CONTROL			###


# This script automates a basic quality control workflow for sequencing data.

# Create a directory to store the quality control reports from FastQC.
mkdir qc-reports

# Run FastQC on all files ending with '.fastq.gz' in the 'data' directory.
# The '-o' flag specifies the output directory for the reports.
fastqc data/*fastq.gz -o qc-reports

# List the newly created reports to verify the command worked.
ls -lh qc-reports

# Note: The next line is for manually copying the HTML reports from a server
# using the 'scp' command. It's for reference and is not part of the script's execution.
# scp username@ipaddress:path/to/your/file/-html ./

# Create a new directory for the consolidated MultiQC report.
mkdir multiqc

# Use MultiQC to aggregate all the individual FastQC reports into a single file.
# The 'qc-reports/' argument tells MultiQC where to find the reports.
# The '-o' flag sets the output directory for the final report.
multiqc qc-reports/ -o multiqc

# List the files in the 'multiqc' directory to confirm the final report was created.
ls -lh multiqc

# Note: This command is for manually opening the final report in a browser.
# You'll need to replace 'path/to/your/file/.html' with the correct path.
# xdg-open path/to/your/file/.html


#build a directory for trimmed files
mkdir trimmed


# Define an array named 'SAMPLES' to store the names of our samples.
# Each sample name will be processed one at a time by the loop.
SAMPLES=(
        "child"
        "father"
)

# Start a 'for' loop to iterate over each item in the SAMPLES array.
for SAMPLE in ${SAMPLES[@]}; do


# Run the fastp command for the current sample.
# We use the 'SAMPLE' variable to construct the input and output file paths.
# The -i and -I flags specify the forward and reverse input FASTQ files.
# The -o and -O flags specify the forward and reverse output files.
# The 'done' keyword marks the end of the loop block.
# Specify html and json files.

fastp \
 -i "data/${SAMPLE}_1.fastq.gz" \
 -I "data/${SAMPLE}_2.fastq.gz" \
 -o "trimmed/${SAMPLE}_1.fastq.gz" \
 -O "trimmed/${SAMPLE}_2.fastq.gz" \
 --html "trimmed/${SAMPLE}.html" \
 --json "trimmed/${SAMPLE}.json"
done

# The script will continue to this point after the loop has processed all samples.
echo "finito, noice"


###				ALIGNMENT			###

SAMPLES=(
        "child"
        "father"
)

# If your reference file dont have the index files with it, you have to index it firstly. I have it, and indexing the grCH38 takes a long time, so I am skipping this step.
# bwa index ref/Homo_sapiens_assembly38.fasta

# Create the output directories if they do not already exist.
mkdir repaired
mkdir mapped


# Start a 'for' loop to iterate through each sample name in the SAMPLES array.
# The 'SAMPLE' variable will hold the current sample name for each iteration.
for SAMPLE in "${SAMPLES[@]}"; do
    # Define the read group information for the current sample.
    # This must be done inside the loop so the 'ID' and 'SM' tags are unique for each sample.
    # GATK requires read group information to process BAM files.
    # The format is @RG followed by tab-separated key-value pairs.
    # ID: unique identifier
    # SM: sample name
    # PL: platform (e.g., ILLUMINA)
    # LB: library name
    READ_GROUP="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:Lib1"

    # Step 1: Repair the paired-end reads.
    # This command uses repair.sh (likely from the BBMap suite) to fix any pairing issues.
    # The output is two paired files and one file for any singleton reads.
    repair.sh in1="trimmed/${SAMPLE}_1.fastq.gz" in2="trimmed/${SAMPLE}_2.fastq.gz" out1="repaired/${SAMPLE}_1.fastq.gz" out2="repaired/${SAMPLE}_2.fastq.gz" outsingle="repaired/${SAMPLE}.single.fq"
        echo "reparing done, moving to mapping"
    # Step 2: Align the repaired reads to the reference genome using BWA-MEM.
    # -t 12: specifies 12 threads for parallel processing, speeding up alignment.
    # -R "$READ_GROUP": adds the read group information directly to the BAM file.
    # The output of 'bwa mem' (SAM format) is piped (|) directly to 'samtools view'.
    bwa mem -t 12 -R "$READ_GROUP" \
        ref/Homo_sapiens_assembly38.fasta \
        "repaired/${SAMPLE}_1.fastq.gz" "repaired/${SAMPLE}_2.fastq.gz" \
    | samtools view -b > "mapped/${SAMPLE}.bam"
done

        echo "finito"


###		MARK DUPLICATES			###

# Create output directories for sorted and marked BAM files.
mkdir sorted
mkdir marked

# Define an array with the names of the samples to be processed.
# This makes it easy to loop through each sample.
SAMPLES=(
        "child"
        "father"
)

for SAMPLE in ${SAMPLES[@]}; do

    # Step 1: Sort the BAM file by coordinate.
    # gatk SortSam is a GATK tool that reorders reads based on their
    # alignment position on the reference genome. This is a crucial
    # pre-processing step required by many GATK tools.
    # -I: specifies the input BAM file.
    # -O: specifies the output BAM file.
    # -SORT_ORDER coordinate: tells the tool to sort the reads by their
    # chromosomal location.

gatk SortSam -I mapped/${SAMPLE}.bam -O sorted/${SAMPLE}-sorted.bam -SORT_ORDER coordinate

    # Provide a message to track the script's progress.
echo "sorting is done"

    # Step 2: Mark duplicate reads.
    # gatk MarkDuplicates identifies and flags reads that are likely
    # PCR duplicates. PCR duplicates can lead to false-positive variant calls,
    # so marking them is an essential step in a GATK pipeline.
    # -I: specifies the input sorted BAM file.
    # -O: specifies the output BAM file with duplicates marked.
    # -M: specifies a file to store a metrics report about the duplicates found.

gatk MarkDuplicates -I sorted/${SAMPLE}-sorted.bam -O marked/${SAMPLE}-marked.bam -M marked/${SAMPLE}-metrics.txt


    # Provide a message to track the script's progress.
echo "marking done"

    # Step 3: Build the BAM index file.
    # gatk BuildBamIndex creates a .bai index file for the final BAM file.
    # An index file allows GATK and other tools to quickly access reads
    # at specific locations in the BAM file without having to read the entire file.
    # This is a required step for GATK's variant calling tools.
gatk BuildBamIndex -I marked/${SAMPLE}-marked.bam -O marked/${SAMPLE}-marked.bai

    # Optional: You can also use samtools to create the index.
    # samtools index marked/${SAMPLE}-marked.bam

done

echo "finito"



### 		BASE QUALITY SCORE RECALIBRATION 		###

# Create a directory to store the output of the Base Quality Score Recalibration (BQSR) steps.
# The 'mkdir' command ensures the directory exists before the script attempts to write files to it.
mkdir bqsr

# Define an array with the names of the samples to be processed.
# This makes it easy to loop through multiple samples without repeating code.
SAMPLES=(

        "child"

        "father"

)

# Begin a 'for' loop to iterate through each sample in the SAMPLES array.
# The 'SAMPLE' variable will hold the name of the current sample.
for SAMPLE in ${SAMPLES[@]}; do

    # Step 1: Run BaseRecalibrator to generate a recalibration table.
    # This tool analyzes patterns of mismatches in the input BAM file against known variant sites.
    # It identifies systematic errors in the base quality scores and creates a model to correct them.
    # -I: Input BAM file, which should be sorted and have duplicates marked.
    # -R: Reference genome.
    # --known-sites: These are crucial. GATK uses these files of known variants (from databases like dbSNP)
    #                to build an accurate recalibration model.
    # -O: Output file for the recalibration table.
gatk BaseRecalibrator \
 -I marked/${SAMPLE}-marked.bam \
 -R ref/Homo_sapiens_assembly38.fasta \
 --known-sites ref/Homo_sapiens_assembly38.dbsnp138.vcf \
 --known-sites ref/Homo_sapiens_assembly38.known_indels.vcf.gz \
 -O bqsr/${SAMPLE}-recal.table

echo "recal table is made, moving to bqsr"

    # Step 2: Apply the recalibration model to create a new BAM file with corrected base quality scores.
    # This step uses the table created by BaseRecalibrator to adjust the quality scores of each base in the reads.
    # -I: Input BAM file (the same one used in the previous step).
    # -R: Reference genome.
    # --bqsr-recal-file: The recalibration table generated in the previous step.
    # -O: The final output BAM file with recalibrated quality scores.
gatk ApplyBQSR \
 -I marked/${SAMPLE}-marked.bam \
 -R ref/Homo_sapiens_assembly38.fasta \
 --bqsr-recal-file bqsr/${SAMPLE}-recal.table \
 -O bqsr/${SAMPLE}-bqsr.bam

done

echo "finito"


### 			VARIANT CALLING				###
# Create a directory to store the output VCF files.
# The 'mkdir' command ensures the directory exists before the script attempts to write files.
mkdir vcf

# Define an array with the names of the samples to be processed.
SAMPLES=(
    "child"
    "father"
)

# Start a 'for' loop to iterate through each sample in the SAMPLES array.
# The 'SAMPLE' variable will hold the name of the current sample for each loop.
for SAMPLE in "${SAMPLES[@]}"; do
    # Run the GATK HaplotypeCaller tool for individual sample variant calling.
    # HaplotypeCaller uses a graph-based approach to call variants with high accuracy.
    # It analyzes the aligned reads in the BAM file to identify and call variants.
    # -I: Input BAM file. This should be the BQSR-recalibrated BAM file for the best results.
    # -R: Reference genome.
    # -O: Output file for the variants. The '.g.vcf.gz' extension indicates a gVCF file.
    # -ERC GVCF: This is a crucial flag. It tells HaplotypeCaller to run in "GVCF mode".
    #            Instead of just reporting variants, it outputs a gVCF file which includes
    #            information about confident reference sites as well. This is essential for
    #            joint genotyping downstream.

    gatk HaplotypeCaller \
        -I bqsr/${SAMPLE}-bqsr.bam \
        -R ref/Homo_sapiens_assembly38.fasta \
        -O vcf/${SAMPLE}.g.vcf.gz \
        -ERC GVCF
done

# The script provides a final message to indicate completion.
echo "Variant calling is done, finito"

