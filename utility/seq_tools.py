import argparse
import subprocess


def raxml_bootstrapped_tree(args):
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -# 20 -s " + args.alignment + " -n T13 "], shell=True)
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -b 12345 -# " + args.bootstraps + " -s " + args.alignment + " -n T14 "], shell=True)
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

if __name__ == '__main__':
    main()