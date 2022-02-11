#!/bin/bash
#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
#SBATCH --mem-per-cpu=2G      # in megabytes, unless unit explicitly stated
#SBATCH --time=20:00:00
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=USERNAME@INSTITUTIONAL_ADDRESS # email address used for event notification
#SBATCH --mail-type=BEGIN,END,FAIL # email on job start, end, and/or failure

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
module load samtools/1.10
module load bedtools2-2.27.1-gcc-8.3.1-g76ahj7
module load deeptools/3.3.0-spack

## Useful shortcuts
export refdir=/your/reference/directory
export workingdir=/your/working/directory

## List of sequences
list=("sample_input1" "sample_inputn" "sample_ip1" "sample_ipn" "sample_neg1" "sample_negn" "control_input1" "control_inputn" "control_ip1" "control_ipn" "control_neg1" "control_negn")


## The commands you want to run

# Create new directories to store files
mkdir beds
mkdir wigs

# map forward and reverse reads to genome

for i in ${list[@]}

do
	echo ${i}
	
	# Create an index of the reference genome in FASTA format so that bedtools can access it quickly
	samtools faidx $refdir/SPECIES_GENOME.INDEX.fa
	
	# Bin map sequences to create histograms
	bedtools genomecov -ibam $workingdir/bowtie/${i}.sorted.bam -bg -g $refdir/SPECIES_GENOME.INDEX.fa.fai > $workingdir/beds/${i}.bedgraph
	
	# Generate a coverage track to visualize in IGV
	bamCoverage -b $workingdir/bowtie/${i}.sorted.bam -o $workingdir/wigs/${i}.bw

	# Note: you should use filtered sequences if duplicates are removed
	# Evaluate whether coverage is good before continuing!
done
       

