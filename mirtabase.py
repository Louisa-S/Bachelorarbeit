import csv

with open("ath_MTI_short.csv") as f:
	mirna = f.readlines()

mirnalist = []
for i in range(0, len(mirna)-1):
	m = mirna[i]
	m = m.split(",")
	mirnalist.append(m)
	
print mirnalist[i]
