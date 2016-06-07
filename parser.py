import sys
from Bio.Emboss.Applications import NeedleCommandline
import Bio
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import RNAAlphabet
from Bio import AlignIO

from random import randint


#import xlwt as xlsw

def filterid():
	ids = dict()
	with open("gene_id.txt", "w") as out:
		
		with open(str(sys.argv[1])) as f:
			next(f)
			for line in f:
				dat = line.split(";")
				if dat[4] not in ids:
					ids[dat[4]] = ""
					
					
		for dat in ids:
			out.write(dat+"\n")
				

'''
waterc = NeedleCommandline()
seqa = Seq("ACUG", RNAAlphabet())
seqb = Seq("AUCGUCGCUAGACUGUUU", RNAAlphabet())
waterc.asequence = "asis:ACUG"
waterc.bsequence = "asis:AUCGGGACUG"
waterc.gapopen = 10
waterc.gapextend = 0.5
waterc.outfile = "out.txt"
print waterc
print waterc.outfile

stdout, stderr = waterc()
print stdout+stderr
'''


def idfile():
	strong, weak, other = 0, 0, 0

	with open(str(sys.argv[1])) as f:
		next(f)
		for line in f:
			
			line = line.split(";")
			s = line[3].split("\r")
			s = s[0].split("\n")
			
			if s[0] == "Functional MTI":
				strong += 1
			else:
				if s[0] == "Functional MTI (Weak)":
					weak += 1
				else:
					other += 1

	print strong
	print weak
	print other


def positives():	
	
	with open(str(sys.argv[1])) as f:
		#next(f)
		mirna = []
	
		#for line in f:
		l = f.readline()
		
		
		while len(mirna) < 1001:	
		
			line = f.readline()	
			if line == "":
				continue
			
			if line[0:3] != "hsa":
				continue	
			
			s = line.split(";")
				
			mirna.append(s[2])
				
	
	with open(str(sys.argv[2])) as g:
		#next(g)
		target = []
		
		l = g.readline()
		#for line in g:
		print l
			
		while len(target) < 1001:	
			line = g.readline()
			
			if line == "":
				continue
			
			s = line.split(";")
				
			target.append(s[1])
				
	
	
	with open("nonfunct.csv", "w") as out:
		
		for i in range(len(mirna)):
			out.write(str(mirna[i])+";"+str(target[i])+";"+"X"+";"+"Non-Functional MTI"+"\n")
		
		
positives()




def table():

	workbook = xlsw.Workbook("hello.xlsx")
	worksheet = workbook.add_worksheet()

	worksheet.write('A1', 'Hello world')

	workbook.close()





