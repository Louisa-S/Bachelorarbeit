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

#command: time python mirna.py utr.txt hsa_MTI_id.csv miRNA_short.csv genid_refseq.txt X X X X

#returns a dictionary with gene ID as key and a list of NM numbers as value
#parameter: name of the file, output of ID converter
#tap limited 
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



#returns a dictionary with miRNA name as key, miRNA sequence as value
#only mature miRNAs are stored
#parameter: name of file from mirBase
#csv file with ; as seperator
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
					
	return mirnalist
	
	
	
#returns a dictionary with miRNA name as value, list of target gene IDs as value
#parameter: name of file from miRTarBase
#csv file with ; as seperator
#excludes anything but Functional MTI
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
	
	
	return mirnalist	
	



#returns a dictionary with NM number as key, list of 5' UTR, gen sequence and 3' UTR as value
#paramter: name of file from UCSC Site
#simple text document, description line starts with >
#UTR lower case, gen upper case
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
			if line == "\n":
				break
				
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
	
	


#for every entry in miRTarBase find match 
#match: Sequence of miRNA, sequence of 5' UTR, gene, 3' UTR 
#paramter: dictionaries of targets, miRNAs, genes and the conversion list	
def findmatch(mirtarbase, mirbase, genes, gentable, p1, p2, p3, p4):
	
	score = []
	pair = []
		
	#matrix with one miRNA per line
	#1 for matching nucleotides in alignment, O for mismatch
	matrix = []
	
	for mlist in mirtarbase:
		
		for mirtar in mirtarbase[mlist]:
			match = []
			
			if mirtar in gentable:
				nm = gentable[mirtar]
			elif mirtar == "":
				continue
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
												
					match.append(genes[n][0])			
					match.append(genes[n][1])			
					match.append(genes[n][2])	
					
						
					#reverse complement of miRNA 
					 
					seqmirna = Seq(match[0], RNAAlphabet())
					seqmirna = seqmirna.reverse_complement()
					
					seq5utr = Seq(match[1]) 
					seqgen = Seq(match[2])
					seq3utr = Seq(match[3])
					
					size5utr = len(match[1])					
					sizegen = len(match[2])					
					size3utr = len(match[3])
					
					#transcribe ro mRNA
					
					seq5utrtransc = seq5utr.transcribe()
					seqgentransc = seqgen.transcribe()
					seq3utrtransc = seq3utr.transcribe()
					
		
					completegen = seq5utrtransc + seqgentransc + seq3utrtransc
					completegen = completegen.upper()
									
					
					#alignment of miRNA and target gene sequence
					#return of alignment: list of sequences(2), score, start, end position
					alignment = pairwise2.align.localms(completegen, seqmirna, p1, p2, p3, p4)
			
			
					#computes a line for the matrix, one line is one startposition of alignment
					
					for parts in alignanalysis(alignment, len(seqmirna), score, sizegen, size5utr, n, mlist):						
						matrix.append(parts)
						
		
	
	
	with open("scores+positions.txt", "w") as f:
		for i in score:
			f.write(str(i)+"\n")	
				
	analysematrix(matrix)				



#computes list with 1 for match and O for mismatch in aligment
#parameter: computed alignment of find match function, size of respective miRNA, list to store scores, size of gene, size of 5' UTR, NM number, miRNA name
def alignanalysis(alignment, mirnasize, score, sizegen, size5, nm, mlist):
	
	startpos = []	
	tupel = ""
	maxi = 22
	
	#list of miRNA name, NM number and alignment score
	tupel += str(mlist)+ " " + str(nm) + " " + str(alignment[0][2]) + " "
	
	#multiple alignments possible because of same score
	for ali in alignment:		
		i = 0
		matrix = []
		
		#computes the start position of alignment
		while (ali[1][i] == "-"):
			i = i+1
		
		if i not in startpos:
			startpos.append(i)
			mirnaseq = ali[1][i:i+mirnasize]
			mseq = ali[0][i:i+mirnasize]
			
			#adds the start position depending on 3' UTR
			#will be negative if not in 3' UTR
			tupel += str(ali[3] - sizegen - size5 + 1)
			tupel += " "
			
			
			#reverse sequences to start from 5' end  
			mirnaseq = mirnaseq[::-1]
			mseq = mseq[::-1]
			
			#computes array of 1 and 0 for the nucleotide positions
			for index in range(min(len(mirnaseq), len(mseq))):
				if mirnaseq[index] == mseq[index]:
					matrix.append(1)
				else:
					matrix.append(0)
						
						
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
	plt.ylim([0,0.8])
	plt.xlabel("Nucleotide position")
	plt.ylabel("Ratio of complementary nucleotides")
	plt.savefig("resultsneu/non-ratio5-1-8-4_new.png")
	
	#stores the percentages for each position in a file
	with open("resultsneu/non-result5-1-8-4_new.txt", "w") as result:
		result.write("Positions"+"\t"+"Ratio of complementary nucleotides"+"\n")
		for i in range(22):
			result.write(str(i)+"\t"+str(perces[i])+"\n")
	
	
	
	
	
	
	
	 
