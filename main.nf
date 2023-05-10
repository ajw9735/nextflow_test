
#! /usr/bin/env nextflow

samples = ["NA18757", "NA18627", "NA18591"]


// Define the input channel for raw fastq files
Channel.fromPath("${params.input_dir}/*_R1.fastq.gz")
    .pair(Channel.fromPath("${params.input_dir}/*_R2.fastq.gz"))
    .filter { it.baseName.startsWith(samples) }
    .set { input_ch }

// Task to run fastp on the raw fastq files
process fastp {
    input:
    file read1 from input_ch
    file read2 from input_ch.pairPath
    
    output:
    file("${output_dir}/${read1.baseName}_processed_R1.fastq.gz") into processed_ch
    file("${output_dir}/${read2.baseName}_processed_R2.fastq.gz") into processed_ch.pairPath
    
    script:
    """
    fastp -i ${read1} -I ${read2} -o ${output}/${read1.baseName}_processed_R1.fastq.gz -O ${output}/${read2.baseName}_processed_R2.fastq.gz
    """
}

// Task to align the processed fastq files with bwa
process bwa_align {
    input:
    file read1 from processed_ch
    file read2 from processed_ch.pairPath
    
    output:
    file("${output_dir}/${read1.baseName}.sam") into aligned_ch
    
    script:
    """
    bwa mem -t 8 ${params.reference_genome} ${read1} ${read2} > ${output}/${read1.baseName}.sam
    """
}

// Task to convert SAM alignments to BAM format with samtools
process sam_to_bam {
    input:
    file sam from aligned_ch
    
    output:
    file("${output_dir}/${sam.baseName}.bam") into bam_ch
    
    script:
    """
    samtools view -bS ${sam} > ${output}/${sam.baseName}.bam
    """
}

// Task to coordinate sort the alignments with Picard SortSam
process sort_bam {
    input:
    file bam from bam_ch
    
    output:
    file("${output_dir}/${bam.baseName}_sorted.bam")
    
    script:
    """
    picard SortSam I=${bam} O=${output}/${bam.baseName}_sorted.bam SORT_ORDER=coordinate
    """
}

// Define the output channels for each task
processed_ch = Channel.create()
aligned_ch = Channel.create()
bam_ch = Channel.create()

// Run the pipeline
fastp | bwa_align | sam_to_bam | sort_bam

// Define the final output channels for the sorted BAM files
output_bams = Channel.from(samples) { "${output_dir}/${it}/${it}_processed_R1_sorted.bam" }

// Copy the final output files to the appropriate output channels
sort_bam.out.into { output_bams[it.baseName.split("_")[0]] }
