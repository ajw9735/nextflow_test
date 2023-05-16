**INSTRUCTIONS:**

This is designed to work using the slurm executor on an HPC. It can also be run locally.
You MUST import a relevant reference genome and name it "Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa".

1. Clone the repo and enter the directory:

`cd nf_fastp_bwa`

2. From the directory import nextflow:

`module load nextflow/21.10.6`

3. Then run this command to see options:

`nextflow run script.nf --help`

4. Then you can run this command to execute the script, changing the username:

`sbatch nextflow run script.nf --usrname name@nyu.edu`

5. Results will show up for each process in the results-dir directory, with the final 3 bam files in the sort_results subdirectory.

