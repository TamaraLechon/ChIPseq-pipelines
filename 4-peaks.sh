#!/bin/bash
#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
#SBATCH --mem-per-cpu=2G      # in megabytes, unless unit explicitly stated
#SBATCH --time=40:00:00
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
module load py-numpy-1.19.0-gcc-8.3.1-gxcejiv
module load py-macs2-2.2.4-gcc-8.3.1-jswbbot


## Useful shortcuts
export refdir=/your/reference/directory
export workingdir=/your/working/directory


## The commands you want to run

# Create a directory to store files
mkdir peaks

# Create a statistical model to define what mapped reads form a peak

# Note1: macs2 -g is effective or mappable genome size, it will vary depending on species of interest
# Note2: macs2 --keep-dup: default behaviour is to only leave one optical duplicate. If picard is run, then --keep-dup should be kept as all
# Note3: macs2 --broad: add flag if analysing histone peaks
# Note4: macs2 -B will create bedGraph files for fragment pileup if building lambda model



# Call peaks for sample ips using corresponding input as control

macs2 callpeak -t $workingdir/bowtie/sample_ip1.sorted.bam -c $workingdir/bowtie/sample_input1.sorted.bam -f BAMPE  \
	       	-g 119146348 --outdir peaks -n sample-ip-rep1 -B

# Repeat n times for each replicate 
macs2 callpeak -t $workingdir/bowtie/sample_ipn.sorted.bam -c $workingdir/bowtie/sample_inputn.sorted.bam -f BAMPE  \
	       	-g 119146348 --outdir peaks -n sample-ip-repn -B


# Call peaks for sample negative controls using corresponding input as control

macs2 callpeak -t $workingdir/bowtie/sample_neg1.sorted.bam -c $workingdir/bowtie/sample_input1.sorted.bam -f BAMPE \
			-g 119146348 --outdir peaks -n sample_neg-rep1 -B

# Repeat n times for each replicate
macs2 callpeak -t $workingdir/bowtie/sample_negn.sorted.bam -c $workingdir/bowtie/sample_inputn.sorted.bam -f BAMPE \
			-g 119146348 --outdir peaks -n sample_neg-repn -B


# Call peaks for control ips or additional sample group/treatment using corresponding input

macs2 callpeak -t $workingdir/bowtie/control_ip1.sorted.bam -c $workingdir/bowtie/control_input1.sorted.bam -f BAMPE  \
	       	-g 119146348 --outdir peaks -n control-ip-rep1 -B

# Repeat n times for each replicate 
macs2 callpeak -t $workingdir/bowtie/control_ipn.sorted.bam -c $workingdir/bowtie/control_inputn.sorted.bam -f BAMPE  \
	       	-g 119146348 --outdir peaks -n control-ip-repn -B


# Call peaks for corresponding negative controls

macs2 callpeak -t $workingdir/bowtie/control_neg1.sorted.bam -c $workingdir/bowtie/control_input1.sorted.bam -f BAMPE \
			-g 119146348 --outdir peaks -n control_neg-rep1 -B

# Repeat n times for each replicate
macs2 callpeak -t $workingdir/bowtie/control_negn.sorted.bam -c $workingdir/bowtie/control_inputn.sorted.bam -f BAMPE \
			-g 119146348 --outdir peaks -n control_neg-repn -B
