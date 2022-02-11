#!/bin/bash
#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=4G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=USERNAME@INSTITUTIONAL_ADDRESS	# email
#SBATCH --mail-type=BEGIN,END,FAIL	# email on job start, end, and/or failure

echo "Usable Environment Variables:"
echo "============================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

## Load some Modules
module load fastqc/v0.11.9
module load fastp/v0.20
module load multiqc/1.9

## Set up working directory
export workingdir=/your/working/dir


## The commands you want to run

# List of sequences
list=("sample_input1" "sample_inputn" "sample_ip1" "sample_ipn" "sample_neg1" "sample_negn" "control_input1" "control_inputn" "control_ip1" "control_ipn" "control_neg1" "control_negn")

# fastqc the raw data (assuming PE data)
for i in ${list[@]}
do
        echo ${i}

	fastqc $workingdir/${i}_R1.fastq.gz
	fastqc $workingdir/${i}_R2.fastq.gz

done

# Summarize QC data of raw sequences
multiqc -i "PROJECT_NAME_RAW_SEQUENCES" $workingdir/

# Trim low quality reads and remove adaptors
for i in ${list[@]}
do
        echo ${i}

	fastp -i $workingdir/${i}_R1.fastq.gz -I $workingdir/${i}_R2.fastq.gz -o $workingdir/trimmed/${i}_fp1.fastq.gz -O $workingdir/trimmed/${i}_fp2.fastq.gz --detect_adapter_for_pe --trim_poly_g --correction

done

# fastqc trimmed sequences
for i in ${list[@]}
do
        echo ${i}

	fastqc $workingdir/trimmed/${i}_fp1.fastq.gz
	fastqc $workingdir/trimmed/${i}_fp2.fastq.gz

done


# Summarize QC data of trimmed sequences
multiqc -i "PROJECT_NAME" $workingdir/trimmed/