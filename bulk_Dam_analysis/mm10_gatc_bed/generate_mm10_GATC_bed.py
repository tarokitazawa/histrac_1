from Bio import SeqIO
import sys

fasta_file = "/path/to/your/BSgenome.Mmusculus.UCSC.mm10.fa"
recognition_sequence = "GATC"  # Modify this as per your enzyme

with open("mm10_GATC.bed", "w") as output:
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq)
        start_pos = seq_str.find(recognition_sequence)
        while start_pos != -1:
            end_pos = start_pos + len(recognition_sequence)
            output.write(f"{record.id}\t{start_pos}\t{end_pos}\n")
            start_pos = seq_str.find(recognition_sequence, start_pos + 1)
