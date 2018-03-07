import os

def fasta_dict(fasta):
    if not os.path.exists(fasta):
        raise Exception('File not found; check path')
    d = {}
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                C_seq = ''
                C_split_line = line.split(' ')
                C_name = C_split_line[0]
                C_name = C_name.rstrip()
                C_name = C_name.lstrip('>')
            else:
                C_seq = C_seq + str(line)
                C_seq = C_seq.rstrip()
                C_seq = C_seq.upper()

            d[C_name] = C_seq
    return(d)