import os
import subprocess

class Utility():

    def fasta_dict(self, fasta):
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

    def rc(self, seq):

        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
        useq = seq.upper()
        letters = list(useq)
        completters = [basecomplement[base] for base in letters]
        rc_list =  completters[::-1]
        rc_seq = "".join(rc_list)
        return(rc_seq)

    def codons(self, seq):
        ORF_list = []
        for i in range(len(seq)):
            codon = [seq[j:j+3] for j in range(i, len(seq), 3)]

            current_seq	 = []
            for i in codon:
                if i != 'TAG' and i != 'TAA' and i != 'TGA':
                    current_seq.append(i)
                else:
                    joined_seq = ''.join(current_seq)
                    ORF_list.append(joined_seq)
                    del(current_seq[:])
                    break

        return(max(ORF_list, key = len))

    def get_longest_orf(self, fasta):
        seqs = self.fasta_dict(fasta)
        stop_codons = ('TAG','TAA','TGA')
        for s in seqs:
            seq = seqs[s]
            compseq = self.rc(s)
            # if no stop codons in the sequence, print the sequence
            if any(stop not in s for stop in stop_codons) or any(stop not in compseq for stop in stop_codons):
                print('This sequence has no stop codons - cannot determine correct reading frame')
                print('>' + s + '\n' + seq + '\n')
            else:
                FandR = []
                F = self.codons(s)
                R = self.codons(self.rc(s))
                FandR.append(F)
                FandR.append(R)
                print('>' + s + '\n' + max(FandR, key = len))



    def reciprocal_best_hits(self):
        pass
        #WRITE


