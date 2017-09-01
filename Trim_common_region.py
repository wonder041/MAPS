from Bio import SeqIO
import multiprocessing
import os.path
import math
import io
import sys

ref_nucl_path, ref_prot_path, input_faa_with_ref_path, input_fna_path, output_faa_path, output_fna_path = sys.argv[1:]

aln_nuc_len = len(next(SeqIO.parse(ref_nucl_path, "fasta")))
PR_nuc_min = aln_nuc_len
PL_nuc_max = 0

for seqr in SeqIO.parse(ref_nucl_path, "fasta"):
    if "_L" in seqr.id:
        PL_nuc = len(str(seqr.seq).rstrip("-"))
        if PL_nuc > PL_nuc_max:
            PL_nuc_max = PL_nuc
    if "_R" in seqr.id:
        PR_nuc = aln_nuc_len - len(str(seqr.seq).lstrip("-"))
        if PR_nuc < PR_nuc_min:
            PR_nuc_min = PR_nuc

tit_PL_PR_tup_dic = {}
for seqr in SeqIO.parse(ref_nucl_path, "fasta"):
    if ("MIMI" in seqr.id) or ("gene" in seqr.id):
        PL = math.ceil(len(seqr.seq[:PL_nuc_max].ungap("-")) / 3)
        PR = len(seqr.seq[:PR_nuc_min].ungap("-")) // 3  # inner+1
        tit_PL_PR_tup_dic[seqr.id] = (PL, PR)


def Gen_pos_in_aln(seq):
    for i in range(len(seq)):
        if seq[i] == "-":
            continue
        yield i

        
        
aln_len = len(next(SeqIO.parse(input_faa_with_ref_path, "fasta")).seq)
PR_min = aln_len
PL_max = 0

for seqr in SeqIO.parse(input_faa_with_ref_path, "fasta"):
    if not seqr.id.startswith("MEGA"):
        break
    pos_in_aln_arr = list(Gen_pos_in_aln(str(seqr.seq)))

    PL = pos_in_aln_arr[tit_PL_PR_tup_dic[seqr.id][0]]
    PR = pos_in_aln_arr[tit_PL_PR_tup_dic[seqr.id][1]]
    if PL > PL_max:
        PL_max = PL
    if PR < PR_min:
        PR_min = PR

# print(PL_max,PR_min)
# sys.exit(0)
###############################################################


def Ge_trim_fna(id):
    skip_l = int(id.split("L")[-1].split("R")[0])
    skip_r = int(id.split("R")[-1])
    frame_shift = int(id.split("F")[-1].split("L")[0].split("C")[0])
    trim_fna = (id.split("F")[0], (frame_shift +
                                   skip_l * 3, frame_shift + skip_r * 3))
    return trim_fna


def Trim(seqr):
    skip_l = len(seqr.seq[:PL_max].ungap("-"))
    skip_r = len(seqr.seq[:PR_min].ungap("-"))
    seqr.seq = seqr.seq[PL_max:PR_min].ungap("-")
    if not len(seqr):
        return None
    seqr.description = seqr.name = seqr.id = seqr.id + \
        'L' + str(skip_l) + 'R' + str(skip_r)
    return(seqr)


def Split(input_faa_with_ref_path, size=10000000):
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
    pos_p, pos = pos_tup
    input_file = open(input_faa_with_ref_path, "r")
    input_file.seek(pos_p)
    handle = io.StringIO(input_file.read(pos - pos_p))
    output_handle = io.StringIO()
    trimed_seqr_gen = (Trim(seqr) for seqr in SeqIO.parse(
        handle, "fasta") if "C" not in seqr.id if not seqr.id.startswith("MIMI") if "gene" not in seqr.id)
    trimed_seqr_arr =[ seqr for seqr in trimed_seqr_gen if seqr]
    trim_fna_dict_seg = [Ge_trim_fna(seqr.id) for seqr in trimed_seqr_arr]
    SeqIO.write(trimed_seqr_arr, output_handle, "fasta")
    return (output_handle.getvalue(), trim_fna_dict_seg)


fileseg_gen = Split(input_faa_with_ref_path)
pool = multiprocessing.Pool(4)
str_gen = pool.imap_unordered(Gen_str, fileseg_gen)

trim_fna_dic = {}
output_file = open(output_faa_path, "w")
for (output_str, trim_fna_dict_seg) in str_gen:
    output_file.write(output_str)
    trim_fna_dic.update(trim_fna_dict_seg)

# print(len(trim_fna_dic))
# 950867
##########################################################################


def Trim_fna(seqr):
    try:
        pos_start, pos_end = trim_fna_dic[seqr.id]
        seqr.seq = seqr.seq[pos_start:pos_end]
        seqr.description = seqr.name = seqr.id = seqr.id + \
            'L' + str(pos_start) + 'R' + str(pos_end)
        return(seqr)
    except KeyError:
        return


fna_seqr_gen = (Trim_fna(seqr)
                for seqr in SeqIO.parse(input_fna_path, "fasta"))
fna_seqr_noNone_gen = (seqr for seqr in fna_seqr_gen if seqr)
SeqIO.write(fna_seqr_noNone_gen, output_fna_path, "fasta")
