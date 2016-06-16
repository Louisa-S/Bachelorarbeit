import sys
from Bio.Emboss.Applications import NeedleCommandline
import Bio
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import RNAAlphabet
from Bio import AlignIO

from random import randint
import matplotlib.pyplot as plt


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
		
		
#positives()




def table():

	workbook = xlsw.Workbook("hello.xlsx")
	worksheet = workbook.add_worksheet()

	worksheet.write('A1', 'Hello world')

	workbook.close()
	

#first argument score+positions, second argument mirtarget_positions
def jointable():
	#check if there are positions in the table
	
	pos = dict()
	
	with open(str(sys.argv[2])) as f:
		for line in f:
			if line == "":
				break
			
			##line = line.replace("\n", " ")
			line = line[:-2].split(" ")
			
			pos[str(line[0])+" "+str(line[1])] = line[2:]
			
	
	scores = []
	
	with open(str(sys.argv[1])) as g:
		
		for line in g:
			if line == "":
				break
					
			line = line[1:-2].split(",")
			
			scores.append(line)
		
	found = False 
	count = 0
	countall = 0
	with open("result_positions.txt", "w") as out:
		
		for s in scores:
			countall = len(scores)

			posis = []
			 
			key = s[0][1:-1] + " " + s[1][2:-1]
			
			out.write(str(s)+" ")
			
			if key in pos:
				posis = pos[key]			
				
				if len(posis) == 0:
					continue				
				else:
					for i in range(3, len(s)):
						if "-" in s[i]:
							continue
						for p in posis:							
							if (int(s[i][1:]) >= int(p)-5) and (int(s[i][1:]) <= int(p)+5):
								count += 1
							
			out.write(str(count)+" / "+str(countall)+" \n")		
			found = False	
		
			
	
jointable()



