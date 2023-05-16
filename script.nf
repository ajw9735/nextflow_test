#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
log.info """
        ===========================================
         F A S T P / B W A    P I P E L I N E
  
         Usage:
        -------------------------------------------
         --outdir             : directory for storage of output variables, default is results-dir
         --samples            : List of sample names, default are those specified in project outline 
         --reference          : Path to reference genome, default is just using it from local directory
         --usrname            : HPC username
        ===========================================

        Example command:
	sbatch nextflow run script.nf --usrname "name@nyu.edu" --outdir results
        """
         .stripIndent()

}

//Define channel
reads_ch = Channel.fromFilePairs(["/scratch/work/courses/BI7653/hw2.2023/${params.samples[0]}_{1,2}.filt.fastq.gz","/scratch/work/courses/BI7653/hw2.2023/${params.samples[1]}_{1,2}.filt.fastq.gz","/scratch/work/courses/BI7653/hw2.2023/${params.samples[2]}_{1,2}.filt.fastq.gz"])

// Define process that calls fastp
process fastp {

    publishDir params.outdir, pattern: 'fastp_trim/*'

    // Input and output declarations
    input:
    tuple val(pair_id), path(reads)

    output:
    file('fastp_trim/trim_*')

    // Command to run fastp
    script:
    """
    #!/bin/bash
    #
    #SBATCH --nodes=4
    #SBATCH --tasks-per-node=8
    #SBATCH --cpus-per-task=8
    #SBATCH --time=12:00:00
    #SBATCH --mem=36GB
    #SBATCH --job-name=fastp
    #SBATCH --mail-type=FAIL
    #SBATCH --mail-user=params.usrname

    module load fastp/intel/0.20.1

    mkdir fastp_trim
    fastp -i ${reads[0]} -I ${reads[1]} -o fastp_trim/trim_${reads[0]} -O fastp_trim/trim_${reads[1]}
    """
}

//process fo alignment and sam generation
process bwa {
    cpus 8
    memory '48GB'
    executor 'slurm'    

    publishDir params.outdir, pattern: 'bwa_results/*'

    // Input and output declarations
    input:
    file(ref)
    file(reads)

    output:
    file('bwa_results/output_*')

    // Command to run bwa
    script:
    """
    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --cpus-per-task=8
    #SBATCH --time=8:00:00
    #SBATCH --mem=64GB
    #SBATCH --job-name=bwa_mem
    #SBATCH --mail-type=FAIL
    #SBATCH --mail-user=params.usrname

    module load bwa/intel/0.7.17    
    module load samtools/intel/1.14
    
    mkdir bwa_results
    cp /scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.amb .
    cp /scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.ann .
    cp /scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.bwt .
    cp /scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.pac .
    cp /scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.sa .


    bwa mem -M -t 1 ${ref} ${reads[0]} ${reads[1]} > ${reads[0]}.sam
    samtools view -Sb -h -f 0x2 ${reads[0]}.sam > bwa_out_${reads[0]}.bam 
    samtools sort bwa_out_${reads[0]}.bam -o bwa_results/output_${reads[0]}.bam
    """
}

process sort {
    cpus 8
    memory '48GB'
    executor 'slurm'

    publishDir params.outdir, pattern: 'sort_results/*'
    //Define inputs
    input:
    file(bam_file)

    output:
    file('sort_results/final_*')
    // Script to run picard tools
    script:
    """
    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --cpus-per-task=8
    #SBATCH --time=8:00:00
    #SBATCH --mem=48GB
    #SBATCH --job-name=fastp
    #SBATCH --mail-type=FAIL
    #SBATCH --mail-user=params.usrname    

    module purge
    
    mkdir sort_results

    module load picard/2.17.11

    java -Xmx42g -XX:ParallelGCThreads=1 -jar "/scratch/ajw9735/nf_fastp_bwa/picard.jar" SortSam \
    -INPUT ${bam_file} \
    -OUTPUT sort_results/final_${bam_file} \
    -SORT_ORDER coordinate \
    -MAX_RECORDS_IN_RAM 10000000 \
    -VALIDATION_STRINGENCY LENIENT
    """
}

// Define the workflow
workflow {
    fastp(reads_ch)
    ref = file(params.reference)
    bwa(ref, fastp.out)
    sort(bwa.out)
}

