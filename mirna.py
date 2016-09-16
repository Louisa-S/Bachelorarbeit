import sys
import analysis
import numpy as np	


#mirtarbase parsen
#mirna db parsen
#utr file parsen	

p1 = sys.argv[6]
p2 = sys.argv[7]
p3 = sys.argv[8]
p4 = sys.argv[9]

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

analysis.findmatch(mirtars, mirnas, genelist, convertedgenes, p1, p2, p3, p4)
