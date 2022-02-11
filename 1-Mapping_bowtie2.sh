#!/bin/bash
#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=8      #   
#SBATCH --mem-per-cpu=4G     # in megabytes, unless unit explicitly stated
#SBATCH --time=48:00:00
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=USERNAME@INSTITUTIONAL_ADDRESS  # email
#SBATCH --mail-type=BEGIN,END,FAIL      # email on job start, end, and/or failure


## Load some modules
module load bowtie2/v2.4.1
module load samtools/1.10
module load bamtools/v2.5.1


## Point to directory containing the reference genome where your sequences will map
export refdir=/your/reference/directory

## Declare your working directory
export workingdir=/your/working/directory


## The commands you want to run

# Index genome for quicker access by bowtie2 during alignment
bowtie2-build $refdir/SPECIES_GENOME.dna.toplevel.fa $refdir/SPECIES_GENOME.INDEX.fa
#Example for A. thaliana below:
#bowtie2-build $refdir/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $refdir/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

# List of sequences to align
list=("sample_input1" "sample_inputn" "sample_ip1" "sample_ipn" "sample_neg1" "sample_negn" "control_input1" "control_inputn" "control_ip1" "control_ipn" "control_neg1" "control_negn")

# Create a new directory to store alignment files
mkdir bowtie

# Map forward and reverse reads to genome

for i in ${list[@]}

do
	echo ${i}
	
	# Align sequences to genome
	bowtie2 --maxins 500 --fr -p 8 -x $refdir/SPECIES_GENOME.INDEX.fa -1 $workingdir/trimmed/${i}_fp1.fastq.gz -2 $workingdir/trimmed/${i}_fp2.fastq.gz -S $workingdir/bowtie/${i}.sam
	
	# Compress sam (alignment files) to bam files (binary files)
	samtools view -bS $workingdir/bowtie/${i}.sam > $workingdir/bowtie/${i}.bam
	
	# Organise mapped reads and index for faster access for downstream processing
	samtools sort -@ ${SLURM_CPUS_PER_TASK} -o $workingdir/bowtie/${i}.sorted.bam $workingdir/bowtie/${i}.bam

	samtools index $workingdir/bowtie/${i}.sorted.bam
	
	# Run some stats on mapping sequences
	bamtools stats -in $workingdir/bowtie/${i}.sorted.bam > $workingdir/bowtie/${i}.sorted.stats.txt
	
	# Note: actually have a look at the mapping stats. If not enough sequences have been mapped, it means the data is too low quality, so any downstream processing will be questionable
	
done
