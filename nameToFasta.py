

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

path = input('give me a file with a list of names ')
names = (path)
path2 = input('give me a fasta file from which to get the sequences ')
fasta = (path2)

namesList = []
f = open(names)
lines = f.readlines()
for l in lines:
	name = l.rstrip('\n')
	namesList.append(name)

out = open('nameToFastaOut.fasta', 'w')

seqDict = fastaDict(fasta)
for i in namesList:
	out.write('>' + i + '\n' + seqDict[i] + '\n')
