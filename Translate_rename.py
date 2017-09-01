import sys
import copy
from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord

input_path = sys.argv[1]
fna_output_path = sys.argv[2]
faa_output_path = sys.argv[3]

barcode = input_path.split('/')[-1].split('.')[0]

#read clstr file and generate dict of OTU representative seq id and OTU size
with open(input_path + ".clstr") as clstr_file:
        size = 0
        rep_size_dic = {}
        next(clstr_file)
        for line in clstr_file:
            if line.strip().endswith('*'):
                rep = line.split('>')[1].split('.')[0]
            if line.startswith('>Cluster'):
                rep_size_dic[rep] = size
                size = 0
                continue
            size += 1
        rep_size_dic[rep] = size

def rename(number, seqr):
    '''reanme seqr by barcode and OTU size'''
    seqr.id = barcode + "N" + \
        str(number) + "S" + str(rep_size_dic[seqr.id])
    return seqr

def translate(seqr):
    '''translate fna seqr by 3 reading frames'''
    for frameshift in range(3):
        new_seqr = copy.copy(seqr)
        new_seqr.id = new_seqr.id + "F" + str(frameshift)
        
        #trim seq by frameshift
        length_remainder=(len(new_seqr.seq) - frameshift) % 3
        if  length_remainder== 0:
            pre_trans = new_seqr.seq[frameshift:]
        else:
            pre_trans = new_seqr.seq[frameshift:-length_remainder]
        
        #translate
        new_seqr.seq = pre_trans.translate()
        
        # discard sequences with stop codon
        if '*' in new_seqr.seq:
            continue
            
        yield new_seqr

#sort seqr
seqr_sorted_arr=sorted(SeqIO.parse(input_path, "fasta"), key=lambda seqr:seqr.id)
        
        
#rename fna seqr
rename_seqr_arr = [rename(number, seqr) for number,
                   seqr in enumerate(seqr_sorted_arr)]

#write renamed fna to file                   
SeqIO.write(rename_seqr_arr, fna_output_path, "fasta")

#translate renamed fna seqr to faa seqr
translated_seqr_gen = (
    translated_seqr for seqr in rename_seqr_arr for translated_seqr in translate(seqr))
    
#write translated faa seqr to file
SeqIO.write(translated_seqr_gen, faa_output_path, "fasta")
