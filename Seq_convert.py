import sys
from Bio import SeqIO

input_path,output_path=sys.argv[1:]

def Detect_format(path):
    extent=path.split(".")[-1]
    if extent in ("faa","fna","fasta","fa"):
        return "fasta"
    if extent in ("fastq",):
        return "fastq"
    if extent in  ("phylip",):
        return "phylip"
    if extent in ("sto","stockholm"):
        return "stockholm"

seqr_gen=SeqIO.parse(input_path,Detect_format(input_path))
SeqIO.write(seqr_gen,output_path,Detect_format(output_path))
