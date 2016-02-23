import csv

with open("miRNA_short.csv") as csvfile:
	mirna = csvfile.readlines()
	
i=0
#list contains a list with 4 elements
#id, precurser, 5p, 3p
mirnalist = []
while i < len(mirna):
	mirna[i] = mirna[i][:-16]
	m = mirna[i].split(",")
	mirnalist.append(m)
	i += 1
	
print mirnalist[1]

