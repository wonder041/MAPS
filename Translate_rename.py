import sys
import copy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Check_OTU_size

input_path = sys.argv[1]
fna_output_path = sys.argv[2]
faa_output_path = sys.argv[3]

barcode = input_path.split('/')[-1].split('.')[0]
rep_size_dic = Check_OTU_size.Get_dic(input_path + ".clstr")

def rename(number, seqr):
    seqr.id = barcode + "N" + \
        str(number) + "S" + str(rep_size_dic[seqr.id])
    return seqr


def translate(seqr):
    for phase in range(3):
        new_seqr = copy.copy(seqr)
        new_seqr.id = new_seqr.id + "F" + str(phase)
        pre_trans = new_seqr.seq[phase:-((len(new_seqr.seq) - phase) % 3)] if (
            len(new_seqr.seq) - phase) % 3 != 0 else new_seqr.seq[phase:]
        trans = pre_trans.translate()
        if '*' in trans:
            # discard sequences with stop codon
            continue
            # new_seqr.id += 'C'
        new_seqr.description = new_seqr.name = new_seqr.id
        new_seqr.seq = trans.ungap('*')
        yield new_seqr

rename_seqr_arr = [rename(number, seqr) for number,
                   seqr in enumerate(SeqIO.parse(input_path, "fasta"))]
SeqIO.write(rename_seqr_arr, fna_output_path, "fasta")
translated_seqr_gen = (
    translated_seqr for seqr in rename_seqr_arr for translated_seqr in translate(seqr))
SeqIO.write(translated_seqr_gen, faa_output_path, "fasta")
