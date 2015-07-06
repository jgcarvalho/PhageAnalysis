from sys import argv

f = open(argv[1])
lines = f.readlines()
pep_i = 1
for i in lines:
	i = i.split()
	if len(i) > 1:
		if len(i[0]) == 16:
			if '*' not in i[0]:
				print ">pep{} [Freq: {}]".format(pep_i, i[1])
				print i[0][2:13]
				pep_i += 1
