import os
import argparse
import pandas
from utility.util import Utility



class Blast():
    def __init__(self, query, db, outfmt=None, max_target_seqs=None, evalue=None, out='out.blast'):
        self.query = query
        self.db = db
        if not outfmt:
            self.outfmt = '6 qseqid sseqid qstart qend sstart send length pident evalue'
        else:
            self.outfmt = outfmt
        self.max_target_seqs = max_target_seqs
        self.evalue = evalue
        self.out = out
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
        self.outfmt = '"' + self.outfmt + '"'
        command = command + ' -outfmt ' + self.outfmt
        if self.max_target_seqs:
            command = command + ' -max_target_seqs ' + str(self.max_target_seqs)
        if self.evalue:
            command = command + ' -evalue ' + str(self.evalue)
        if self.out:
            command = command + ' -out ' + str(self.out)
        return command

    def blastn(self):
        return self.blast_command('blastn')

    def blastp(self):
        return self.blast_command('blastp')

    def blastx(self):
        return self.blast_command('blastx')

    def tblastn(self):
        return self.blast_command('tblastn')

    def db_hits_fasta(self):
        out_fields = self.outfmt.split()
        out_fields.pop(0)
        sseq = out_fields.index('sseqid')
        csv = pandas.read_csv(self.out, sep='\t', header=None)
        db_names = list(csv[sseq])
        return db_names

    def write_hits_fasta(self):
        u = Utility()
        db_seqs = u.fasta_dict(self.db)
        db_names = self.db_hits_fasta
        fasta_out = self.out + '.fasta'
        with open('fasta_out') as outfile:
            for n in db_names:
                outfile.write('>' + n + '\n' + db_seqs[n] + '\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', required=True)
    parser.add_argument('--db', required=True)
    parser.add_argument('--program', required=True)
    parser.add_argument('--outfmt', required=False)
    parser.add_argument('--max_target_seqs', required=False)
    parser.add_argument('--evalue', required=False)
    parser.add_argument('--out', required=False)
    args = parser.parse_args()

if __name__ == '__main__':
    main()
