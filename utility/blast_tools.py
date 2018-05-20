import os
import argparse
import pandas
import subprocess
from utility.util import Utility
from utility.seq_tools import SeqTools


class Blast():
    def __init__(self, query, db, outfmt, out, max_target_seqs, evalue):
        self.outfmt = outfmt
        self.evalue = evalue
        self.out = out
        self.query = query
        self.db = db
        self.max_target_seqs = max_target_seqs
        self.check_files()

    @property
    def check_db(self):
        prot_db_suffix = ['.phr', '.pin', '.psq']
        nucl_db_suffix = ['.nhr', '.nin', '.nsq']
        prot_db_check = list(set([os.path.exists(self.db + i) for i in prot_db_suffix]))[0]
        nucl_db_check = list(set([os.path.exists(self.db + i) for i in nucl_db_suffix]))[0]
        if os.path.isfile(self.db) and prot_db_check or nucl_db_check:
            return 0
        else:
            print('Your database %s either does not exist or it has not been formatted - check full path and db format' % self.db)
            return 1

    @property
    def check_query(self):
        if os.path.isfile(self.query):
            return 0
        else:
            print('Your query file %s is not found - did you enter the full path?' % self.query)
            return 1

    def check_files(self):
        return self.check_db + self.check_query

    def blast_command(self, program):
        command = '%s -query %s -db %s' % (program, self.query, self.db)
        command = command + ' -outfmt ' + self.outfmt
        if self.max_target_seqs:
            command = command + ' -max_target_seqs ' + str(self.max_target_seqs)
        if self.evalue:
            command = command + ' -evalue ' + str(self.evalue)
        if self.out:
            command = command + ' -out ' + str(self.out)
        return command

    def mafft_command(self):
        command = 'mafft %s > %s' % (self.out + '.fasta', self.out + '.fasta.aln')
        return command

    def db_hits_fasta(self):
        out_fields = self.outfmt.split()
        out_fields.pop(0)
        sseq = out_fields.index('sseqid')
        csv = pandas.read_csv(self.out, sep='\t', header=None)
        db_names = list(set(list(csv[sseq])))
        return db_names


    def write_hits_fasta(self, args):
        min_length = args.min_hit_length or 0
        u = Utility()
        db_seqs = u.fasta_dict(self.db)
        db_names = self.db_hits_fasta()
        fasta_out = self.out + '.fasta'
        with open(fasta_out, 'w') as outfile:
            for n in db_names:
                if len(db_seqs[n]) > min_length:
                    outfile.write('>' + n + '\n' + db_seqs[n] + '\n')

    def align(self):
        subprocess.call([self.mafft_command()], shell=True)

    def run(self, args):
        command = self.blast_command(args.program)
        subprocess.call([command], shell=True)
        with open(self.out) as outfile:
            num_hits = outfile.readlines()
            if len(num_hits) > 0:
                self.write_hits_fasta(args)
                if len(num_hits) > 1:
                    self.align()
                else:
                    print('Not enough hits for an alignment')
            else:
                print('No hits found')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', required=True)
    parser.add_argument('--db', required=True)
    parser.add_argument('--program', required=True)
    default_outfmt = '"6 qseqid sseqid qstart qend sstart send length pident evalue"'
    parser.add_argument('--outfmt', default=default_outfmt)
    parser.add_argument('--max_target_seqs', default=30)
    parser.add_argument('--evalue', default='1e-5')
    parser.add_argument('--out', default='out.blast')
    parser.add_argument('--min_hit_length', default=0)
    args = parser.parse_args()
    b = Blast(args.query, args.db, args.outfmt, args.out, args.max_target_seqs, args.evalue)
    b.run(args)

if __name__ == '__main__':
    main()
