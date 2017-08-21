def Decide_mode(read_path):
    with open(read_path, "r") as read_file:
        next(read_file)
        line = next(read_file)
        if "PP" in line:
            return True
        else:
            return False


def Get_dic(read_path):
    if Decide_mode(read_path):
        Measure_line = lambda line: int(line.split('S')[1].split('F')[0])
    else:
        Measure_line = lambda line: 1

    with open(read_path) as read_file:
        size = 0
        tit_size_dic = {}

        next(read_file)
        for line in read_file:
            if '*' in line:
                tit = line.split('>')[1].split('.')[0]
            if line.startswith('>'):
                tit_size_dic[tit] = size
                size = 0
                continue
            size += Measure_line(line)
        tit_size_dic[tit] = size
        return(tit_size_dic)

if __name__ == "__main__":
    import sys
    tit_size_dic = Get_dic(sys.argv[1])
    for tit in list(tit_size_dic.keys())[:10]:
        print(tit, tit_size_dic[tit])
