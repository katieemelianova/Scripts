import subprocess

def percentFloat(proportion,total):
	percentage = (float(proportion) * 100)/float(total)
	return(percentage)

def fastaDict(i):
	f = open(i)
	L = {}

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

		L[C_name] = C_seq
	return(L)


path = input('give me the path to you transcriptome file ')
t = (path)

path2 = input('give me a path to an A. thaliana trancriptome file ')
a = (path2)

transcriptome = open(t)
transcriptomeDict = fastaDict(trancriptome)