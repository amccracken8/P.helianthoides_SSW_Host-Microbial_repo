#!/bin/sh

# Specify a partition
#SBATCH --partition=bluemoon
# Request nodes
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8  # 8 threads per task
# Request memory
#SBATCH --mem=50G
# Run for some time
#SBATCH --time=24:00:00 
# Name job:
#SBATCH --job-name=star_quant2
# Name output file
#SBATCH --output=%x_%j.out

# change to the directory where you submitted this script
cd ${SLURM_SUBMIT_DIR}

# Executable section: echoing some Slurm data
echo "Starting sbatch script myscript.sh at:`date`"


#path to code to be run
cd home/data/AK_pycno_RNA/star

# run code

for i in $(ls home/data/AK_pycno_RNA/01_QC_FastP/fastp_cleaned |cut -f 1 -d "_"| uniq);
do

    echo "starting sample ${i}"

    read1=home/data/AK_pycno_RNA/01_QC_FastP/fastp_cleaned/${i}_1_clean.fq.gz
    read2=home/data/AK_pycno_RNA/01_QC_FastP/fastp_cleaned/${i}_2_clean.fq.gz
    gtfreference=home/data/AK_pycno_RNA/pycno_ref/pycno.gtf
	
    STAR --genomeDir home/data/AK_pycno_RNA/star/pycno_genome_index/ --runThreadN 8 --readFilesCommand gunzip -c --readFilesIn ${read1} ${read2} --outFileNamePrefix home/data/AK_pycno_RNA/star/results_quant_2/${i}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile $gtfreference --sjdbGTFtagExonParentTranscript Parent

    echo "sample ${i} done"
	
done


