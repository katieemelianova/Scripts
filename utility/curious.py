import argparse

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    tenx = subparsers.add_parser('10x')
    tenx.add_argument('--genome_size', required=True, help='haploid genome size in Gbp')
    tenx.set_defaults(func=tenx_spex)
    args = parser.parse_args()
    args.func(args)

def tenx_spex(args):
    bp = (float(args.genome_size) * 1000000)
    num_reads = (bp * 56)/150
    string = 'For a haploid genome size of %s, ' \
             'you will need to generate %s reads to ' \
             'achieve 56x raw covergae' % (args.genome_size, int(num_reads))
    print(string)

if __name__ == '__main__':
    main()
