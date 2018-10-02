from Bio import SeqIO
import multiprocessing
import os.path
import math
import io
import sys

ref_nucl_path, input_faa_with_ref_path, input_fna_path, output_faa_path, output_fna_path = sys.argv[1:]

#calculate alignment length and set initial value of left and right boundary in nucleotide file
aln_nuc_len = len(next(SeqIO.parse(ref_nucl_path, "fasta")))
right_boundary_nuc = aln_nuc_len
left_boundary_nuc = 0

#update boundary from primer positions
for seqr in SeqIO.parse(ref_nucl_path, "fasta"):
    if "_L" in seqr.id:
        left_pos_nuc = len(str(seqr.seq).rstrip("-"))
        if left_pos_nuc > left_boundary_nuc:
            left_boundary_nuc = left_pos_nuc
    if "_R" in seqr.id:
        right_pos_nuc = aln_nuc_len - len(str(seqr.seq).lstrip("-"))
        if right_pos_nuc < right_boundary_nuc:
            right_boundary_nuc = right_pos_nuc

#calculate left and right boundary in non-alignment position for all reference
tit_boundary_dic = {}
for seqr in SeqIO.parse(ref_nucl_path, "fasta"):
    if seqr.id.startswith("MEGA"):
        left_char_pos = math.ceil(len(seqr.seq[:left_boundary_nuc].ungap("-")) / 3)
        right_char_pos = len(seqr.seq[:right_boundary_nuc].ungap("-")) // 3  # inner+1
        tit_boundary_dic[seqr.id] = (left_char_pos, right_char_pos)


def Gen_char_pos_in_aln(seq):
    for i in range(len(seq)):
        if seq[i] == "-":
            continue
        yield i

        
#calculate alignment length and set initial value of left and right boundary        
aln_len = len(next(SeqIO.parse(input_faa_with_ref_path, "fasta")).seq)
right_boundary = aln_len
left_boundary = 0


#left and right boundary in non-alignment position back to alignment position
for seqr in SeqIO.parse(input_faa_with_ref_path, "fasta"):
    if seqr.id.startswith("MEGA"):
        pos_in_aln_arr = list(Gen_char_pos_in_aln(str(seqr.seq)))
    
        left_pos = pos_in_aln_arr[tit_boundary_dic[seqr.id][0]]
        right_pos = pos_in_aln_arr[tit_boundary_dic[seqr.id][1]]
        if left_pos > left_boundary:
            left_boundary = left_pos
        if right_pos < right_boundary:
            right_boundary = right_pos
    else:
        break

        
###############################################################

def Trim(seqr):
    # if seqr.id.startswith("MEGA"):
        # return None
    left_char_pos = len(seqr.seq[:left_boundary].ungap("-"))
    right_char_pos = len(seqr.seq[:right_boundary].ungap("-"))
    seqr.seq = seqr.seq[left_boundary:right_boundary].ungap("-")
    if not len(seqr):
        return None
    if seqr.id.startswith("MEGA"):
        return seqr
    seqr.description = seqr.name = seqr.id = seqr.id + \
        'L' + str(left_char_pos) + 'R' + str(right_char_pos)
    return(seqr)


def Split_input_file(input_faa_with_ref_path, size=10000000):
    input_file_b = open(input_faa_with_ref_path, "rb")
    file_size = os.path.getsize(input_faa_with_ref_path)
    pos_p = 0
    while pos_p + size < file_size:
        input_file_b.seek(pos_p + size)
        while True:
            line = input_file_b.readline()
            if line.startswith(b">"):
                pos = input_file_b.tell() - len(line)
                yield (pos_p, pos)
                pos_p = pos
                break
    yield (pos_p, file_size)


def Gen_str(pos_tup):
    
    #load file segment into handle
    pos_p, pos = pos_tup
    input_file = open(input_faa_with_ref_path, "r")
    input_file.seek(pos_p)
    handle = io.StringIO(input_file.read(pos - pos_p))
    
    #set output to output_handle
    output_handle = io.StringIO()
    
    #trim faa seqr and fliter empty result
    trimed_seqr_gen = (Trim(seqr) for seqr in SeqIO.parse(handle, "fasta"))        
    trimed_seqr_arr =[ seqr for seqr in trimed_seqr_gen if seqr]
    
    #generate dict for trim fna
    # id_pos_tup_arr = [Get_id_pos_tup(seqr.id) for seqr in trimed_seqr_arr]
    
    SeqIO.write(trimed_seqr_arr, output_handle, "fasta")
    return output_handle.getvalue()


file_segment_gen = Split_input_file(input_faa_with_ref_path)
pool = multiprocessing.Pool(4)
str_gen = pool.imap_unordered(Gen_str, file_segment_gen)

trim_fna_dic = {}
output_file = open(output_faa_path, "w")
for output_str in str_gen:
    output_file.write(output_str)
    # for line in output_str.split("\n"):
    #     if line.startswith(">"):
    #         trim_fna_dic[line.split("F")[0].strip(">")]=line.split("F")[1]


##########################################################################
sys.exit(0)

def Trim_fna(seqr):
    try:
        id_ext = trim_fna_dic[seqr.id]

        left_pos = int(id_ext.split("L")[-1].split("R")[0])
        right_pos = int(id_ext.split("R")[-1])
        frame_shift = int(id_ext.split("F")[-1].split("L")[0].split("C")[0])

        seqr.seq = seqr.seq[frame_shift + left_pos * 3: frame_shift + right_pos * 3]
        seqr.description = seqr.name = seqr.id = seqr.id +"F"+ id_ext
        return(seqr)
    except KeyError:
        return


fna_seqr_gen = (Trim_fna(seqr)
                for seqr in SeqIO.parse(input_fna_path, "fasta"))
fna_seqr_noNone_gen = (seqr for seqr in fna_seqr_gen if seqr)
SeqIO.write(fna_seqr_noNone_gen, output_fna_path, "fasta")
