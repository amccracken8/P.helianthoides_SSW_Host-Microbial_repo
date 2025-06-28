from Bio import SeqIO
import re

def subset_fasta(input_fasta, gene_list, output_fasta):
    selected_genes = set()
    
    # Read gene patterns from the list
    with open(gene_list, 'r') as gene_file:
        for line in gene_file:
            gene_pattern = line.strip()
            # Modify the pattern to include 't#' at the end
            gene_pattern += r'\.t\d+$'
            # Match gene IDs with the provided pattern
            selected_genes.update([record.id for record in SeqIO.parse(input_fasta, 'fasta') if re.match(gene_pattern, record.id)])

    # Write selected sequences to the output file
    with open(output_fasta, 'w') as output_file:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.id in selected_genes:
                SeqIO.write(record, output_file, 'fasta')

# Example usage
input_fasta = '/users/a/r/armccrac/data/AK_pycno_RNA/pycno_ref/Pycno_coding_sequences.fasta'
gene_list = 'star_deseqfilt_transcript_ids.txt'
output_fasta = 'star_deseqfilt.transcripts.fasta'

subset_fasta(input_fasta, gene_list, output_fasta)