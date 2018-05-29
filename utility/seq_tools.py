import argparse
import subprocess


def raxml_bootstrapped_tree(args):
    alignment = args.alignment
    num_bootsraps = args.bootstraps
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -# 20 -s " + alignment + " -n T13 "], shell=True)
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -b 12345 -# " + num_bootsraps + " -s " + alignment + " -n T14 "], shell=True)
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15"], shell=True)

def bmge(args):
    default_outfile = args.alignment + ".bmge"
    subprocess.call(["java -jar ~/Desktop/Software/BMGE-1.12/BMGE.jar -i " + args.alignment + " -t " + args.seqtype + " -of " + default_outfile], shell=True)

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    raxml = subparsers.add_parser('raxml')
    raxml.add_argument('--alignment', required=True)
    raxml.add_argument('--bootstraps', required=True)
    raxml.set_defaults(func=raxml_bootstrapped_tree)

    bmge = subparsers.add_parser('bmge')
    bmge.add_argument('--alignment', required=True)
    bmge.add_argument('--seqtype', required=True, help="AA,DNA,RNA,CODON")
    bmge.add_argument('--outfile', required=False)
    bmge.set_defaults(func=bmge)

    args = parser.parse_args()
    args.func(args)

def run_augustus():
    #wrappeer function for a basic augustus run
    pass

def parse_augustus_outfile(annotation_file):

    outfile = annotation_file + '.fasta'
    with open(annotation_file) as openfile:
        line_number = len(openfile.readlines())
    f = open(annotation_file, 'r')
    annotated_seqs = {}
    gene_name = ''
    seq = ''

    for i in range(line_number):
        line = f.next()
        if line.startswith('# start gene'):
            gene_name = line.split()[-1]
        if line.startswith('# coding sequence = '):
            addline = line.lstrip('# ')
            addline = addline.lstrip('coding sequence = [').rstrip('\n')
            seq+=str(addline)
            continue
        if line.startswith('# ') and seq and ']' not in line:
            addline = line.lstrip('# ').rstrip('\n')
            seq+=str(addline)
        elif line.endswith(']\n') and seq:
            addline = line.lstrip('# ').rstrip(']\n')
            seq+=str(addline)
            annotated_seqs[gene_name] = seq
            gene_name = ''
            seq = ''

    with open(outfile, 'w') as o:
        for seq in annotated_seqs:
            o.write('>' + seq + '\n' + annotated_seqs[seq] + '\n')

if __name__ == '__main__':
    main()