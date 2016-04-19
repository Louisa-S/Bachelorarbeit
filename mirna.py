import sys
import analysis


#mirtarbase parsen
#mirna db parsen
#utr file parsen	
		
convertedgenes = analysis.symbolconverter(str(sys.argv[4]))
print convertedgenes[0]

#filter out non human targets
#filter out duplicates with different experiment
#filter out weak ones
mirtars = analysis.mirtarparser(str(sys.argv [2]))
print mirtars[0]
mirnas = analysis.mirbaseparser(str(sys.argv [3]))
print mirnas[0]
genelist = analysis.parsegenes(str(sys.argv [1]))
print genelist[0]

print ("data parsing was successful!")

analysis.findmatch(mirtars, mirnas, genelist, convertedgenes)
