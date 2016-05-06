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

#command: time python mirna.py utr.txt hsa_MTI_short2.csv miRNA_short.csv ucsc_symbols_nm.txt 

#returns list of nm numbers and gene symbols
def symbolconverter(filename):
	table = dict()
	result = dict()
	gene = []
	with open(filename, "r") as f:
		next(f)
		for s in f:
			s1 = s.split("\n")
			s1 = s1[0].split("\r")
			s1 = s1[0].split("\t")
			
			if s1[1] not in table:
				nm = []
				nm.append(s1[0])				
				table[s1[1]] = nm
		
			else:
				nm = []
				nm = table[s1[1]]
				nm.append(s1[0])
				table[s1[1]] = nm 
				
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
				
				mirnalist[m[2]] = m[3]
				
				if m[4] != "":
					mirnalist[m[4]] = m[5]
	
	#maxi = len(mirnalist[0][5])
	#
	#for m1 in mirnalist:
	#	if maxi < len(m1[5]):
	#		maxi = len(m1[5])
	#
	#print maxi		
	return mirnalist

def mirtarparser(filename):
	mirnalist = []
	with open(filename) as f:
		next(f)
		for mirna in f:
			m = mirna.split("\n")
			m = m[0].split("\r")
			m = m[0].split(";")
			
			mirnalist.append(m)
		
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
	
	for m1 in mirbase:
		if maxi < len(mirbase[m1]):
			maxi = len(mirbase[m1])
	
	print maxi	
	
	#matrix with one mirna per line
	#X for match in alignment, O for mismatch
	matrix = []
	for mirtar in mirtarbase:
		match = []
		found = False
		
		#there can be more than 1 NM number 
		if mirtar[1] in gentable:
			nm = gentable[mirtar[1]]
		else:
			print "Gensymbol not found: "+mirtar[1]
			#print mirtar[1]
			continue
			
		#print mirtar[0]
		#print mirtar[0][:-3]
		
		if mirtar[0].lower() in mirbase:
			mi = mirbase[mirtar[0].lower()]
			#print mi
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
			
		for n in nm:
			if n in genes:
				
				match.append(genes[n][0])			
				match.append(genes[n][1])			
				match.append(genes[n][2])	
				
				#print match		
			 
				seqmirna = Seq(match[0], RNAAlphabet())
				seq5utr = Seq(match[1]) 
				size5utr = len(match[1])
				seqgen = Seq(match[2])
				sizegen = len(match[2])
				seq3utr = Seq(match[3])
				size3utr = len(match[3])
				#print sizegen
				seq5utrtransc = seq5utr.transcribe().complement()
				seqgentransc = seqgen.transcribe().complement()
				seq3utrtransc = seq3utr.transcribe().complement()
				
				#align mirna to genes
				#align only seed sequence to utr+gene+utr
				completegen = seq5utrtransc + seqgentransc + seq3utrtransc
				completegen = completegen.upper()
								
				#seed = seqmirna[-8:-1]
				
				#for a in pairwise2.align.localxs(completegen, seqmirna, -10, -0.5):
				#	print(format_alignment(*a))
				#	exit()
				
				#returns alignment: list of Sequences(2), score, start, end position
				alignment = pairwise2.align.localxs(completegen, seqmirna, -10, -0.5)
				#print alignment
				
				#print alignment[0][1]
				
				#computes a line for the matrix, one line is one startposition of alignment
				for parts in alignanalysis(alignment, len(seqmirna), maxi):
					
					matrix.append(parts)
			
		
			
		
	analysematrix(matrix, maxi)				

		
		
		
		
		
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
				#compare alignment position to size of 5utr, gen and 3utr
						
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
	
	
	
	for a in alignment:
		i = 0
		matrix = []
		
		while (a[1][i] == "-"):
			i = i+1
		
		if i not in startpos:
			startpos.append(i)
			mirnaseq = a[1][i:i+mirnasize]
			mseq = a[0][i:i+mirnasize]
		
			#print mirnaseq
			#print mseq
		
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
	
	for x in range(maxlength+1):
		xaxis.append(x)
	
	
	#for every character in the sequence, assuming 28, if mirnalength is < 28 it is filled with gaps "-" at the end
	for i in range(maxlength+1):	
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
	
	plt.show()
	
	
	
	
	
	
	
	 
