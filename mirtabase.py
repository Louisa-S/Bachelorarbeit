import csv

#returns list of tupels with lines of table
def dbparser(filename):
	with open(filename) as f:
		mirna = f.readlines()

	mirnalist = []
	for m in mirna:
		m = m.split(";")
		mirnalist.append(m)
		
	return mirnalist



#filter out gene symbols,
#filter out dubplicates, write list of symbols to uniquesymbols	

def filtergenymbols(filename, mirnas):
	writefile = open("genesymbols.txt", "w")

	for m in mirnas:
		s = m[1]
		writefile.write(s)
		writefile.write("\n")
	
	writefile.close()	


	with open("genesymbols.txt") as sym:
		lines = sym.readlines()
		lines_set = set(lines)
		
		unique = open("uniquesymbols.txt", "w")
		
		for l in lines_set:
			unique.write(l)
			
		unique.close()
