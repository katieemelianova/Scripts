from utility import util, blast_tools

class TestUtil():

    def test_fasta_dict(self):
        pass


class TestBlastTools():
    def __init__(self):
        self.default_outfmt = '"6 qseqid sseqid qstart qend sstart send length pident evalue"'
        self.max_target = 1
        self.e = 1e-5
        self.out = 'out.blast'
        self.db = '/Users/katieemelianova/PycharmProjects/Scripts/tests/test_files/test_db.fasta'
        self.query = '/Users/katieemelianova/PycharmProjects/Scripts/tests/test_files/test_query.fasta'
        self.b = blast_tools.Blast(self.query,
                                   self.db,
                                   outfmt=self.default_outfmt,
                                   max_target_seqs=self.max_target,
                                   evalue=self.e,
                                   out=self.out)


        self.bad_db = '/non/existent/db.fasta'
        self.bad_query = '/non/existent/query.fasta'
        self.bad_b = blast_tools.Blast(self.bad_query,
                                       self.bad_db,
                                       outfmt=self.default_outfmt,
                                       max_target_seqs=self.max_target,
                                       evalue=self.e,
                                       out=self.out)


    def test_check_db(self):
        assert self.b.check_db == 0
        assert self.bad_b.check_db == 1

    def test_check_query(self):
        assert self.b.check_query == 0
        assert self.bad_b.check_query == 1

    def test_check_files(self):
        assert self.b.check_files() == 0
        assert self.bad_b.check_files() == 2

    def test_blast_command(self):
        self.b = blast_tools.Blast(self.query,
                                   self.db, outfmt=self.default_outfmt,
                                   max_target_seqs=self.max_target,
                                   evalue=self.e,
                                   out=self.out)
        assert self.b.blast_command('blastx') == 'blastx ' \
                                                 '-query ' \
                                                 '/Users/katieemelianova/PycharmProjects/Scripts/tests/test_files/test_query.fasta ' \
                                                 '-db /Users/katieemelianova/PycharmProjects/Scripts/tests/test_files/test_db.fasta ' \
                                                 '-outfmt "6 qseqid sseqid qstart qend sstart send length pident evalue" ' \
                                                 '-max_target_seqs 1 ' \
                                                 '-evalue 1e-05 ' \
                                                 '-out out.blast'

        self.b = blast_tools.Blast(self.query,
                                   self.db,
                                   outfmt='"6 qseqid sseqid length pident evalue"',
                                   max_target_seqs=2,
                                   evalue=1e-5,
                                   out='outfile')
        assert self.b.blast_command('blastx') == 'blastx ' \
                                                 '-query /Users/katieemelianova/PycharmProjects/Scripts/tests/test_files/test_query.fasta ' \
                                                 '-db /Users/katieemelianova/PycharmProjects/Scripts/tests/test_files/test_db.fasta ' \
                                                 '-outfmt "6 qseqid sseqid length pident evalue" ' \
                                                 '-max_target_seqs 2 ' \
                                                 '-evalue 1e-05 ' \
                                                 '-out outfile'



    def test_db_hits_fasta(self):
        b = blast_tools.Blast(self.query,
                                   self.db,
                                   outfmt='6 qseqid sseqid length pident evalue',
                                   max_target_seqs=2,
                                   evalue=1e-5,
                                   out='/Users/katieemelianova/PycharmProjects/Scripts/tests/test_files/out.blast')
        f = b.db_hits_fasta()

        assert f == ['AT5G13930.1']




