import csv
import Bio
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import RNAAlphabet
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
from itertools import groupby
from Bio.SubsMat.MatrixInfo import pam60

#command: time python mirna.py utr.txt hsa_MTI_short2.csv miRNA_short.csv ucsc_symbols_nm.txt 

#returns list of nm numbers and gene symbols
def symbolconverter(filename):
	table = dict()
	#result = dict()
	gene = []
	with open(filename, "r") as f:
		next(f)
		for s in f:		
			s = s.replace(" ", "")
			s1 = s.split("\t")

			#s1 = s1[0].split("\r")
			#s1 = s1[0].split("\t")
			
			if s1[0] not in table:
				liste = []
				nm = s1[1].split(";")
				for x in nm:
					if x[0:2] == "NM":
						liste.append(x)
								
				table[s1[0]] = liste
		
			#else:
			#	nm = []
			#	nm = table[s1[0]]
			#	nm.append(s1[0])
			#	table[s1[1]] = nm 
				
	#table2 = sorted(table, key = lambda symbol: symbol[1])
		
	#for key, group in groupby(table2, lambda x: x[1]):
	#	gene.append(key)		
	#	nmlist = []
	#	for nm in group:
	#		nmlist.append(nm[0])
	#		
	#	gene.append(nmlist)		
	#	result.append(gene)
	#	gene = []
	#print table	
	return table



#returns list of tupels with lines of table
def mirbaseparser(filename):
	mirnalist = dict()
	with open(filename) as f:
		next(f)
		for mirna in f:
			if mirna[0:3] == "hsa":
				m = mirna.split("\n")
				m = m[0].split("\r")
				m = m[0].split(";")
				
				mirnalist[m[2].lower()] = m[3]
				
				if m[4] != "":
					mirnalist[m[4].lower()] = m[5]
					
	#maxi = len(mirnalist[0][5])
	#
	#for m1 in mirnalist:
	#	if maxi < len(m1[5]):
	#		maxi = len(m1[5])
	#
	#print maxi		
	return mirnalist

def mirtarparser(filename):
	mirnalist = dict()
	#genesymbols = dict()
	with open(filename) as f:
		next(f)
		for mirna in f:
			#if "(Weak)" in mirna:
			#	continue
			m = mirna.split("\n")
			m = m[0].split("\r")
			m = m[0].split(";")
			
			if m[3] != "Functional MTI":
				continue
			#if m[1] not in genesymbols:
			#	genesymbols[m[1]] = 1
				
			if m[0] in mirnalist:
				if m[1] not in mirnalist[m[0]]:
					genes = mirnalist[m[0]]
					genes.append(m[1])
					mirnalist[m[0]] = genes
			else:
				genes = []
				genes.append(m[1])
				mirnalist[m[0]] = genes
			
			#mirnalist.append(m)
	
	#with open("geneconversion.txt", "w") as f:
	#	for x in genesymbols:
	#		f.write(x)
	#		f.write("\n")
			
	return mirnalist	
	


def parsegenes (filename):
	with open(filename) as gene:
		genelist = dict()
		gen = []
		genexist = False
		utr5 = ""
		utr3 = ""
		genseq = ""
		first = 1
		#utr5, utr3, genseq as strings, boolean genexist if gene is between utrs
		for line in gene:
			if first == 1:
				gen.append(line)
				first += 1
				continue
			else:
				if line[0] == ">":
					gen.append(utr5)
					gen.append(genseq)
					gen.append(utr3)
					
					nm = gen[0].split(" ")
		
					nm = nm[0][14:]
							
					genelist[nm] = gen[1:]	
				
					
					gen = []
					gen.append(line)
					utr5 = ""
					utr3 = ""
					genseq = ""
					genexist = False

					#skip to next line
					continue
				else: 
					for i in range(0, len(line)-1):
						if line[i].islower():
							if genexist:
								utr3 += line[i]
							else:
								utr5 += line[i]
						else:
							genseq += line[i]
							genexist = True	
		
		gen.append(utr5)
		gen.append(genseq)
		gen.append(utr3)
		
		nm = gen[0].split(" ")
		
		nm = nm[0][14:]
				
		genelist[nm] = gen[1:]	

	
	return genelist
	

#build reverse complement of genes !!!		
def findmatch(mirtarbase, mirbase, genes, gentable):
	
	maxi = 0
	notfound = 0
	aligncount = 0
	basecount = 0
	ca, cc, cg, cu = 0, 0, 0, 0
	
	for m1 in mirbase:
		if maxi < len(mirbase[m1]):
			maxi = len(mirbase[m1])
	
	print maxi	
	
	#matrix with one mirna per line
	#X for match in alignment, O for mismatch
	matrix = []
	for mlist in mirtarbase:
		
		for mirtar in mirtarbase[mlist]:
			match = []
			found = False
			
			#there can be more than 1 NM number 
			if mirtar in gentable:
				nm = gentable[mirtar]
			else:
				print "Gensymbol not found: "+mirtar
				#print mirtar[1]
				notfound += 1
				continue
				
			#print mirtar[0][:-3]
			
			if mlist.lower() in mirbase:
				mi = mirbase[mlist.lower()]
				
				for i in range(len(mi)-2):
					basecount += 1
					if mi[i:i+2] == "AA":
						ca += 1
						#basecount += 1
					else:
						if mi[i:i+2] == "CC":
							cc += 1
							#basecount += 1
						else:
							if mi[i:i+2] == "GG":
								cg += 1
								#basecount += 1
							else: 
								if mi[i:i+2] == "UU":
									cu += 1
									#basecount += 1
				
				match.append(mi)
				#if mi[1] == mirtar[0].lower():
				#	match.append(mi[2])
				#else:
				#	if mi[3] != "":
				#		match.append(mi[4])
			else:
				continue
			#else:
			#	if mirtar[0][:-3].lower() in mirbase:
			#		mi = mirbase[mirtar[0][:-3].lower()]
			#		if mi[1] == mirtar[0].lower():
			#			match.append(mi[2])
			#		else:
			#			if mi[3] != "":
			#				match.append(mi[4])
			'''
			for n in nm:
				
				if n in genes:
					
					#print mlist
					#print n
					
					#if mlist != "hsa-miR-203a-3p" or n != "NM_007313":
					#	continue
					#print match[0]
								
					match.append(genes[n][0])			
					match.append(genes[n][1])			
					match.append(genes[n][2])	
					
					#print match		
					#reverse complement of mirna 
					#transcribe gene, all Ts to Us
					#alignment with gene
				 
					seqmirna = Seq(match[0], RNAAlphabet())
					seqmirna = seqmirna.reverse_complement()
					seq5utr = Seq(match[1]) 
					size5utr = len(match[1])
					seqgen = Seq(match[2])
					sizegen = len(match[2])
					seq3utr = Seq(match[3])
					size3utr = len(match[3])
					#print sizegen
					#print size3utr
					seq5utrtransc = seq5utr.transcribe()
					seqgentransc = seqgen.transcribe()
					seq3utrtransc = seq3utr.transcribe()
					
					#align mirna to genes
					#align only seed sequence to utr+gene+utr
					completegen = seq5utrtransc + seqgentransc + seq3utrtransc
					completegen = completegen.upper()
									
									
					#for a in pairwise2.align.localms(completegen,seqmirna,1,-2, -2, -1):
					#	print(format_alignment(*a))
	
					
					#returns alignment: list of Sequences(2), score, start, end position
					alignment = pairwise2.align.localms(completegen, seqmirna,1, -2, -2, -1)
					#print alignment
					
					#print alignment[0][1]
					#aligncount += len(alignment)
					#print len(alignment)
					
					#computes a line for the matrix, one line is one startposition of alignment
					#for parts in alignanalysis(alignment, len(seqmirna), maxi):
						
					#	matrix.append(parts)
					'''
				
	
	acontent = float(ca) / float(basecount)
	ccontent = float(cc) / float(basecount)
	gcontent = float(cg) / float(basecount)
	ucontent = float(cu) / float(basecount)
	
	print acontent
	print ccontent	
	print gcontent	
	print ucontent	
	print basecount
	
	print len(matrix)
	#analysematrix(matrix, maxi)				
	#print notfound
	#print aligncount	
	with open("dinucleo_strong.txt", "w") as out:
		out.write("AA content: "+str(acontent)+"\n")
		out.write("CC content: "+str(ccontent)+"\n")
		out.write("GG content: "+str(gcontent)+"\n")
		out.write("UU content: "+str(ucontent)+"\n")
		
		#out.write(str(aligncount))
	#print "number of alignments: "
	
		
		
		
		
	"""
	
	for mi in mirbase:
		if mi[2] in mirtar[0].lower():
			#print mirtar[0]
			#print mi[2]
			#print "found two matching mirna names"
			match.append(mi[3])
			found = True
		else:
			if mi[4] != "":
				if mi[4] in mirtar[0].lower():
					#print mirtar[0]
					#print mi[4]
					#print "found two matching mirna names"
					match.append(mi[5])
					found = True
				# nm = converted nm number
		
		if found:	
			for n in nm:
				if len(match) > 1:
					match.pop()
					match.pop()
					match.pop()
				#print nm
				#print gen[0]
				# if first line starting with > contains nm number
				for gen in genes:							
					if n in gen[0]:#.find(nm) != -1:	#result list or just one?!?		
						match.append(gen[1]) #utr5p
						match.append(gen[2]) #exon
						match.append(gen[3]) #utr3p
						
						#print "found matching refseq number"
						#print gen[0]
						#print nm
						
						# further align this match, reverse complement of genes!! 
						#print match
						
					
						#build complements of genes	
						
							
							#print matrix
							
						
						#analysematrix(matrix)				
						
						#exit()
						
						break
			break
			
			#break"""
		
				#consider experiment type
				
"""	

def genconverter(symbol, gentable):
	for gen in gentable:
		if symbol == gen[0]:
			#print gen[1] 
			#print "converted"
			return gen[1]
	
	return "Gensymbol not found!"
"""


#computes list with X for match and O for mismatch in aligment
def alignanalysis(alignment, mirnasize, maxi):
	#print "aligner"
	
	startpos = []	
	
	
	
	a = alignment[0]
	
	i = 0
	matrix = []
	
	while (a[1][i] == "-"):
		i = i+1
	
	if i not in startpos:
		startpos.append(i)
		mirnaseq = a[1][i:i+mirnasize]
		mseq = a[0][i:i+mirnasize]
	
		#print mirnaseq
		#print mseq.back_transcribe()
	
		for index in range(mirnasize):
			if mirnaseq[index] == mseq[index]:
				matrix.append("X")
			else:
				matrix.append("O")
					
		
		matrix.reverse()

		
		gaps = mirnasize
		while gaps < maxi+1:
			matrix.append("-")
			gaps = gaps + 1
		#print matrix

		yield matrix
	
	
	
	
	
	
def analysematrix(matrix, maxlength):
		
	#the dimension is the number of mirnas analysed
	
								
	perces = []
	xaxis = []
	
	for x in range(maxlength):
		xaxis.append(x)
	
	
	#for every character in the sequence, if mirnalength is < maxlength of a mirna it is filled with gaps "-" at the end
	for i in range(maxlength):	
		sumx = 0 	
		dimension = len(matrix)
		
		for line in matrix:
			if line[i] == "X":
				sumx = sumx + 1
			else:
				if line[i] == "-":
					dimension = dimension - 1
			
		if dimension > 0:
			perc = float(sumx) / float(dimension)
	
			#print perc	
		
			perces.append(perc)
		else:
			perces.append(0)
		
		
	plt.bar(xaxis, perces)
	plt.xlabel("positions in miRNA")
	plt.ylabel("percentage of complementary bases")
	plt.savefig("probability_positions_strong.png")
	
	with open("results_strong.txt", "w") as result:
		result.write("Positions"+"\t"+"Percentage of complementary bases"+"\n")
		for i in range(maxlength):
			result.write(str(i)+"\t"+str(perces[i])+"\n")
	
	
	
	
	
	
	
	 
