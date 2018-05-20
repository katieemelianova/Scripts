import argparse
import subprocess


def raxml_bootstrapped_tree(args):
    alignment = args.alignment
    num_bootsraps = args.bootstraps
    print("~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -# 20 -s " + alignment + " -n T13 ")
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -# 20 -s " + alignment + " -n T13 "], shell=True)
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -b 12345 -# " + num_bootsraps + " -s " + alignment + " -n T14 "], shell=True)
    subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15"], shell=True)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    raxml = subparsers.add_parser('raxml')
    raxml.add_argument('--alignment', required=True)
    raxml.add_argument('--bootstraps', required=True)
    raxml.set_defaults(func=raxml_bootstrapped_tree)

    args = parser.parse_args()
    args.func(args)



if __name__ == '__main__':
    main()