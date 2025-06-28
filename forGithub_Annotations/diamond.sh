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
#SBATCH --job-name=star_all_diamond
# Name output file
#SBATCH --output=%x_%j.out

# change to the directory where you submitted this script
cd ${SLURM_SUBMIT_DIR}

# Executable section: echoing some Slurm data
echo "Starting sbatch script myscript.sh at:`date`"


#path to code to be run
cd /users/a/r/armccrac/data/AK_pycno_RNA/04_diamond

#run code

# for DEGs 
diamond blastx -d /users/a/r/armccrac/scratch/nr/nr_diamond_db.dmnd -q star_deseqfilt.transcripts.fasta --outfmt 6 qseqid sseqid staxids pident length evalue bitscore stitle --max-target-seqs 1 --evalue 0.05 -o NE_Star_ALL_dmnd.tsv

#for echinoderms
#diamond blastx -d /users/a/r/armccrac/Databases/nr/nr_diamond_db.dmnd -q /users/a/r/armccrac/data/AK_pycno_RNA/03_CDHit/fp_cdhit.fasta --taxonlist 7586 --outfmt 6 qseqid sseqid staxids pident length evalue bitscore --max-target-seqs 5 --evalue 0.05 -o echino_diamond_out.tsv

#for bacteria
#diamond blastx -d /users/a/r/armccrac/Databases/nr/nr_diamond_db.dmnd -q /users/a/r/armccrac/data/AK_pycno_RNA/03_CDHit/fp_cdhit.fasta --taxon-exclude 2759 --outfmt 6 qseqid sseqid staxids pident length evalue bitscore --max-target-seqs 5 --evalue 0.05 -o prok_diamond_out.tsv