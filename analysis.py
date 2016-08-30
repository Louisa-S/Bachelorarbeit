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

#command: time python mirna.py utr.txt hsa_MTI_id.csv miRNA_short.csv genid_refseq.txt

#returns a dictionary with gene id as key and a list of nm numbers as values
def symbolconverter(filename):
	table = dict()
	gene = []
	with open(filename, "r") as f:
		next(f)
		for s in f:	
			if s == "\n":
				break
					
			s = s.replace(" ", "")
			s1 = s.split("\t")
				
			if s1[0] not in table:
				liste = []
				nm = s1[1].split(";")
				for x in nm:
					if x[0:2] == "NM":
						liste.append(x)
								
				table[s1[0]] = liste
		
	return table



#returns a dictionary with mirna name as key, mirna sequence as value
#only mature mirnas are stored
def mirbaseparser(filename):
	mirnalist = dict()
	with open(filename) as f:
		next(f)
		for mirna in f:
			if mirna == "\n":
				break
				
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
	
	
#returns a dictionary with minra name as value, list of gene ids as value
def mirtarparser(filename):
	mirnalist = dict()
	with open(filename) as f:
		next(f)
		for mirna in f:
			if mirna == "\n":
				break
				
			m = mirna.split("\n")
			m = m[0].split("\r")
			m = m[0].split(";")
			
			if m[3] != "Functional MTI":
				continue
				
			if m[0] in mirnalist:
				if m[1] not in mirnalist[m[0]]:
					genes = mirnalist[m[0]]
					genes.append(m[1])
					mirnalist[m[0]] = genes
			else:
				genes = []
				genes.append(m[1])
				mirnalist[m[0]] = genes
				
	#with open("geneconversion.txt", "w") as f:
	#	for x in genesymbols:
	#		f.write(x)
	#		f.write("\n")
	
	return mirnalist	
	


#returns a dictionary with NM number as key, list of 5p utr, gen and 3p utr as value
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
			#if line == "\n":
			#	break
				
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
	
	

#for every entry in mirtarbase finds match 
#match: Sequence of mirna, sequence of 5p utr, gene, 3p utr 	
def findmatch(mirtarbase, mirbase, genes, gentable):
	
	maxi = 0	
	#basecount = 0
	#ca, cc, cg, cu = 0, 0, 0, 0
	score = []
	pair = []
	
	#computes the maximum length of all mirnas
	for m1 in mirbase:
		if maxi < len(mirbase[m1]):
			maxi = len(mirbase[m1])
	
	#matrix with one mirna per line
	#1 for match in alignment, O for mismatch
	matrix = []
	
	for mlist in mirtarbase:
		
		for mirtar in mirtarbase[mlist]:
			match = []
			found = False
			
			if mirtar in gentable:
				nm = gentable[mirtar]
			else:
				print "Gensymbol not found: "+mirtar	
				continue				
			
			if mlist.lower() in mirbase:
				mi = mirbase[mlist.lower()]
								
				match.append(mi)
				
			else:
				continue
			
			origmatch = match
			for n in nm:
				
				match = origmatch
				
				if n in genes:
										
					#if mlist != "hsa-miR-106a-5p" or n != "NM_005734":
					#	continue					
								
					match.append(genes[n][0])			
					match.append(genes[n][1])			
					match.append(genes[n][2])	
					
						
					#pair.append((mlist,n))
			
					#reverse complement of mirna 
					#transcribe gene, all Ts to Us
					#alignment with gene
				 
					seqmirna = Seq(match[0], RNAAlphabet())
					seqmirna = seqmirna.reverse_complement()
					
					seq5utr = Seq(match[1]) 
					seqgen = Seq(match[2])
					seq3utr = Seq(match[3])
					
					size5utr = len(match[1])					
					sizegen = len(match[2])					
					size3utr = len(match[3])
					
					seq5utrtransc = seq5utr.transcribe()
					seqgentransc = seqgen.transcribe()
					seq3utrtransc = seq3utr.transcribe()
					
					#align mirna to genes
					
					completegen = seq5utrtransc + seqgentransc + seq3utrtransc
					completegen = completegen.upper()
									
				
					#returns alignment: list of Sequences(2), score, start, end position
					alignment = pairwise2.align.localms(completegen, seqmirna,5,-4,-6,-4)
			
					
					#print alignment[0][1]
					#aligncount += len(alignment)
					#print len(alignment)
					
					#computes a line for the matrix, one line is one startposition of alignment
					
					for parts in alignanalysis(alignment, len(seqmirna), score, sizegen, size5utr, n, mlist):						
						matrix.append(parts)
						
		
	
	
	with open("resultsneu/non-scores+positions5-4-6-4_new.txt", "w") as f:
		for i in score:
			f.write(str(i)+"\n")	
				
	analysematrix(matrix)				



#computes list with 1 for match and O for mismatch in aligment
def alignanalysis(alignment, mirnasize, score, sizegen, size5, nm, mlist):
	
	startpos = []	
	tupel = ""
	maxi = 22
	
	#list of mirnaname, NM number and alignment score
	tupel += str(mlist)+ " " + str(nm) + " " + str(alignment[0][2]) + " "
	
	for ali in alignment:		
		i = 0
		matrix = []
		
		#computes the startposition of alignment
		while (ali[1][i] == "-"):
			i = i+1
		
		if i not in startpos:
			startpos.append(i)
			mirnaseq = ali[1][i:i+mirnasize]
			mseq = ali[0][i:i+mirnasize]
			
			#adds the starposition depending on 3p utr
			#will be negative if not in 3p utr
			tupel += str(ali[3] - sizegen - size5 + 1)
			tupel += " "
			
			mirnaseq = mirnaseq[::-1]
			mseq = mseq[::-1]
			
			#computes array of 1 and 0 for the base positions
			for index in range(min(len(mirnaseq), len(mseq))):
				if mirnaseq[index] == mseq[index]:
					matrix.append(1)
				else:
					matrix.append(0)
						
			#reverses this array to get it from the 5p end on
			#matrix.reverse()
			
			#adds - to make every mirna as long as the longest			
			gaps = mirnasize
			
			while gaps < maxi:
				matrix.append("-")
				gaps = gaps + 1
					
			yield matrix
			
	#adds the tupel with the information to the score list for the final table
	score.append(tupel)	
	
	
	
#makes the plot of the complement base pairs in the alignment
def analysematrix(matrix):	
							
	perces = []
	xaxis = []
	
	#defines the x axis
	for x in range(22):
		xaxis.append(x)
	
	
	#for every character in the sequence, if mirnalength is < maxlength of a mirna it is filled with gaps "-" at the end
	for i in range(22):	
		sumx = 0 	
		
		#the dimension is the number of mirnas analysed		
		dimension = len(matrix)
		
		for line in matrix:
			if line[i] == 1:
				sumx = sumx + 1
			
				#if line[i] == "-":
				#	dimension = dimension - 1
			
		if dimension > 0:
			perc = float(sumx) / float(dimension)		
			perces.append(perc)
		else:
			perces.append(0)
		
	
	#plots the graph
	plt.bar(xaxis, perces)
	plt.xlabel("Nucleotide position")
	plt.ylabel("Ratio of complementary nucleotides")
	plt.savefig("resultsneu/non-ratio5-4-6-4_new.png")
	
	#stores the percentages for each position in a file
	with open("resultsneu/non-result5-4-6-4_new.txt", "w") as result:
		result.write("Positions"+"\t"+"Ratio of complementary nucleotides"+"\n")
		for i in range(22):
			result.write(str(i)+"\t"+str(perces[i])+"\n")
	
	
	
	
	
	
	
	 
