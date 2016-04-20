import csv
import Bio
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import RNAAlphabet
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align import AlignInfo

from itertools import groupby


#returns list of nm numbers and gene symbols
def symbolconverter(filename):
	table = []
	result = []
	gene = []
	with open(filename, "r") as f:
		next(f)
		for s in f:
			s1 = s.split("\n")
			s1 = s1[0].split("\r")
			s1 = s1[0].split("\t")
			
			table.append(s1)
			
	table2 = sorted(table, key = lambda symbol: symbol[1])
		
	for key, group in groupby(table2, lambda x: x[1]):
		gene.append(key)		
		nmlist = []
		for nm in group:
			nmlist.append(nm[0])
			
		gene.append(nmlist)		
		result.append(gene)
		gene = []
		
	return result



#returns list of tupels with lines of table
def mirbaseparser(filename):
	mirnalist = []
	with open(filename) as f:
		next(f)
		for mirna in f:
			if mirna[0:3] == "hsa":
				m = mirna.split("\n")
				m = m[0].split("\r")
				m = m[0].split(";")
			
				mirnalist.append(m)
		
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
		genelist = []
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
					genelist.append(gen)
					
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
		genelist.append(gen)	

	
	return genelist
	

#build reverse complement of genes !!!		
def findmatch(mirtarbase, mirbase, genes, gentable):
	
	for mirtar in mirtarbase:
		match = []
		found = False
		
		#there can be more than 1 NM number 
		nm = genconverter(mirtar[1], gentable)
		
		if nm == "Gensymbol not found!":
			continue
		
		for mi in mirbase:
			if mi[2] in mirtar[0].lower():
				#print mirtar[0]
				#print mi[2]
				print "found two matching mirna names"
				match.append(mi[3])
				found = True
			else:
				if mi[4] != "":
					if mi[4] in mirtar[0].lower():
						#print mirtar[0]
						#print mi[4]
						print "found two matching mirna names"
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
											
							seed = seqmirna[-8:-1]
							
							#for a in pairwise2.align.localxs(completegen, seqmirna, -10, -0.5):
							#	print(format_alignment(*a))
							#	exit()
							
							#returns alignment: list of Sequences, score, start, end position
							alignment = pairwise2.align.localxs(completegen, seed, -10, -0.5)
							print alignment
							
							print alignment[0][1]
							
							#if alignment[0][3] <= size5utr:
							#	print "5 utr"
							#else:
							#	if alignment[0][3] <= sizegen:
							#		print "in gene"
							#	else:
							#		print "3 utr"
							
							exit()
							#print alignment
							break
				break
				
				#break
			
				#consider experiment type
				#compare alignment position to size of 5utr, gen and 3utr
				
						
						


def genconverter(symbol, gentable):
	for gen in gentable:
		if symbol == gen[0]:
			#print gen[1] 
			#print "converted"
			return gen[1]
	
	return "Gensymbol not found!"


