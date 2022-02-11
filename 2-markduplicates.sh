#!/bin/bash
#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=4      #  
#SBATCH --mem-per-cpu=8G       # in megabytes, unless unit explicitly stated
#SBATCH --time=6:0:0
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


## Load some modules
module load picard/2.22.2
module load bamtools/v2.5.1

## Useful shortcuts
export refdir=/your/reference/directory
export workingdir=/your/working/directory

## List of sequences
list=("sample_input1" "sample_inputn" "sample_ip1" "sample_ipn" "sample_neg1" "sample_negn" "control_input1" "control_inputn" "control_ip1" "control_ipn" "control_neg1" "control_negn")

## The commands you want to run

# Create a directory for storing filtered sequences
mkdir markdup


# Filter sequences to leave only uniquely mapped reads
for i in ${list[@]}

do
        echo ${i}

	##  MARK DUPLICATES  ##
	java -jar $PICARD MarkDuplicates I=$workingdir/bowtie/${i}.sorted.bam O=$workingdir/markdup/${i}.markdup.bam M=$workingdir/markdup/${i}.metrics.markdup.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

	bamtools stats -in $workingdir/markdup/${i}.markdup.bam > $workingdir/markdup/${i}.markdup.dupstats.txt


	## REMOVE DUPLICATES ##
	java -jar $PICARD MarkDuplicates I=$workingdir/bowtie/${i}.sorted.bam O=$workingdir/markdup/${i}.rmdup.bam M=$workingdir/markdup/${i}.metrics.rmdup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

	bamtools stats -in $workingdir/markdup/${i}.rmdup.bam > $workingdir/markdup/${i}.rmdup.dupstats.txt

	# Note: actually look at the stats to decide whether it is better to remove or keep duplicates.
	
done
