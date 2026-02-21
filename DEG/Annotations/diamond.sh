#!/bin/sh

# Specify a partition
#SBATCH --partition=bluemoon
# Request nodes
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1 # tasks per node
# Request memory
#SBATCH --mem=50G
# Run for some time
#SBATCH --time=24:00:00 
# Name job:
#SBATCH --job-name=diamond
# Name output file
#SBATCH --output=%x_%j.out

# change to the directory where you submitted this script
cd ${SLURM_SUBMIT_DIR}

# Executable section: echoing some Slurm data
echo "Starting sbatch script myscript.sh at:`date`"


#path to code to be run
cd Home/data/AK_pycno_RNA/04_diamond

#run code

# for DEGs 
diamond blastx -d Home/scratch/nr/nr_diamond_db.dmnd -q star_deseqfilt.transcripts.fasta --outfmt 6 qseqid sseqid staxids pident length evalue bitscore stitle --max-target-seqs 1 --evalue 0.05 -o NE_Star_ALL_dmnd.tsv
