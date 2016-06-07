import sys
import analysis
import numpy as np	


#mirtarbase parsen
#mirna db parsen
#utr file parsen	

'''
with open("david_converted_txt", "w") as output:

	with open("david_symbols.txt") as f:
		next(f)
		for line in f:
			sym = line.split("\t")
			print sym
			if sym[2] == "Homo sapiens":
				output.write(sym[1]+"\t"+sym[0]+"\n")
			else:
				continue
			
'''


convertedgenes = analysis.symbolconverter(str(sys.argv[4]))
#print convertedgenes
print "Converting table done"

#filter out non human targets
#filter out duplicates with different experiment
#filter out weak ones
#accept files with blank lines at the end!!!!
mirtars = analysis.mirtarparser(str(sys.argv [2]))
#print mirtars
print "Mirtar done"
mirnas = analysis.mirbaseparser(str(sys.argv [3]))
#print mirnas
print "Mirna done"
genelist = analysis.parsegenes(str(sys.argv [1]))
#print genelist
print "Genelist done"

#print ("data parsing was successful!")

analysis.findmatch(mirtars, mirnas, genelist, convertedgenes)
