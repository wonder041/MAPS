import sys

tit_path,input_fna_path,fna_tit_path=sys.argv[1:]

with open(tit_path,"r") as tit_file:
    tit_set=frozenset(line.split("F")[0] for line in tit_file)

with open(input_fna_path,"r") as input_fna_file:
    with open(fna_tit_path,"w") as fna_tit_file:
        for line in input_fna_file:
            if not line.startswith(">"):
                continue
            query_tit=line[1:].split("L")[0]
            if query_tit in tit_set:
                fna_tit_file.write(line[1:])
        