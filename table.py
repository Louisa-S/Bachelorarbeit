import sys
import xlsxwriter
import os
import statistics


#python table.py "scores+positions(1-2-2-1).txt" "scores+positions(2-2-5-4).txt" "scores+positions(3-2-4-4).txt" "scores+positions(5-1-8-4).txt" "scores+positions(5-2-8-3).txt" "scores+positions(5-3-8-2).txt" "scores+positions(5-4-6-4).txt"

def firstcol(w, f):
	row = 1
	col = 0
	
	for line in f:
		line = line.split(" ")		
		
		w.write(row, col, line[0]+" - "+line[1])
		
		row += 1 
	
	w.write(row+1, col, "Average")
	w.write(row+2, col, "Standard deviation")
	w.write(row+4, col, "T-Test")
	w.write(row+6, col, "Negative Control")
	

def firstcolnon(w, f, rowlen):
	
	row = rowlen
	col = 0
	
	for line in f:
		line = line.split(" ")		
		
		w.write(row, col, line[0]+" - "+line[1])
		
		row += 1 
	
	w.write(row+1, col, "Average")
	w.write(row+2, col, "Standard deviation")
	w.write(row+4, col, "T-Test")

def firstline(k, w, name):	

	w.write(0, k, "Alignment Score " +str(k)+": "+ name[16:-4] + " / found in miRTarBase")

			

def writecol(i, f, w):
	print i
	row = 1
	col = i
	
	scores = []
	scos = []
		
	for line in f:
		if line == "":
			break
		
		line = line[:-1].split(" ")
		
		scores.append(line)
		scos.append(float(line[2]))
		
	count = 0
	countall = len(scores)
	
	scos = list(scos)
	
	m = statistics.mean(scos)
	st = statistics.stdev(scos)
		
	w.write(countall+2, col, m)
	w.write(countall+3, col, st)
		
		
	for s in scores:	
		
		posis = []
		endstr = ""
		 
		key = s[0] + " " + s[1]
		
		endstr += str(s[2])+" "		
		
		if key in pos:
			posis = pos[key]	
			
			if len(posis) == 0:
				endstr += "("+str(count)+" / "+str(countall)+")"		
				w.write(row, col, endstr)
				row += 1
				continue				
			else:
				for i in range(3, len(s)):
					if s[i] == "":
						continue
					if "-" in s[i]:
						continue
					for p in posis:							
						if (int(s[i]) >= int(p)-5) and (int(s[i]) <= int(p)+5):
							count += 1
						
		endstr += "("+str(count)+" / "+str(countall)+")"
		
		w.write(row, col, endstr)
		row += 1	



def writecolnon(i, f, w, rowpos):
	print i
	row = rowpos
	col = i
	
	scores = []
	scos = []
	
	for line in f:
		if line == "":
			break
		
		line = line[:-1].split(" ")
		
		scores.append(line)
		scos.append(float(line[2]))
		
	count = 0
	countall = len(scores)
	
	scos = list(scos)
		
	m = statistics.mean(scos)
	st = statistics.stdev(scos)
		
	w.write(countall+1 +rowpos, col, m)
	w.write(countall+2 +rowpos, col, st)		
		
	for s in scores:	
		
		posis = []
		endstr = ""
		 
		key = s[0] + " " + s[1]
		
		endstr += str(s[2])+" "		
			
		w.write(row, col, endstr)
		row += 1	

def ratioline(w, k, name):
	w.write(0, k, "Ratio compl bases: "+name[6:-4])
	
	
def ratios(w, f, i):
	print i
	ratios = []
	rats = []
	
	for line in f:
		if line == "":
			break
		else:
			line = line.split("\t")
			ratios.append(line[1][:-1])
			rats.append(float(line[1][:-1]))
	
	row = 1
	for r in ratios:
		w.write(row, i, r)
		row += 1
	
	w.write(row, i, statistics.mean(rats))
	
	w.write(row+2, 0, "Negative control")

def nonratios(w, f, i, filelen):
	print i
	
	ratios = []
	rats = []
	
	for line in f:
		if line == "":
			break
		else:
			line = line.split("\t")
			ratios.append(line[1][:-1])
			rats.append(float(line[1][:-1]))
	
	row = 1 + filelen
	
	for r in ratios:
		w.write(row, i, r)
		row += 1
	w.write(row, i, statistics.mean(rats))
	
			
#main

workbook = xlsxwriter.Workbook('result_table.xlsx')
workbook2 = xlsxwriter.Workbook('ratio_table.xlsx')

w = workbook.add_worksheet()
w2 = workbook2.add_worksheet()

w.write(0, 0, "MTI")
w2.write(0, 0, "Nucleotides")

c = 0
i = 0
for r in range(1,29):
	w2.write(r, c, i)
	i += 1

w2.write(r+1, 0, "Average")

i = 0	
for r in range(33, 61):
	w2.write(r, c, i)
	i += 1 
	
w2.write(r+1, 0, "Average")

#first argument is path of directory with results
#second one is parsed website with positions

pos = dict()
	
with open(str(sys.argv[2])) as f:
	for line in f:
		if line == "":
			break

		line = line[:-2].split(" ")
		
		pos[str(line[0])+" "+str(line[1])] = line[2:]
		
			
path = str(sys.argv[1])
first = True

filelist = os.listdir(path)


count = 0

for i in range(len(filelist)):
	if filelist[i][0:6] == "scores":
		count += 1
	
k = 1
j = 1
targetlen = 0
ratiolen = 0

for i in range (0, len(filelist)):	
	if first:
		if filelist[i][0:6] == "scores":			
			with open(path+"\\"+filelist[i]) as one:
				targetlen = sum(1 for line in one)	
			with open(path+"\\"+filelist[i]) as one2:
				firstcol(w, one2)
			
			first = False
			
	if filelist[i][0:6] == "scores":	
		with open(path+"\\"+filelist[i]) as f:				
			firstline(k, w, str(filelist[i]))	
			writecol(k, f, w)
			k += 1	
			
	if filelist[i][0:6] == "result":
		with open(path+"\\"+filelist[i]) as f2:
			ratiolen = sum(1 for line in f2)
		with open(path+"\\"+filelist[i]) as f:
			next(f)
			ratioline(w2, j, str(filelist[i]))
			ratios(w2, f, j)
			j += 1
		
			
			
k = 1
j = 1
first = True

nontargetlen = 0
nonratiolen = 0

for i in range(len(filelist)):
	if first:
		if filelist[i][0:10] == "non-scores":			
			with open(path+"\\"+filelist[i]) as f:
				nontargetlen = sum(1 for line in f)	
			with open(path+"\\"+filelist[i]) as g:
				firstcolnon(w, g, targetlen + 9)
			
			first = False
		
	if filelist[i][0:10] == "non-scores":	
		with open(path+"\\"+filelist[i]) as g:	
			writecolnon(k, g, w, targetlen + 9)			
			k += 1		
			
	if filelist[i][0:10] == "non-result":
		print j
		with open(path+"\\"+filelist[i]) as f2:
			nonratiolen = sum(1 for line in f2)
		with open(path+"\\"+filelist[i]) as f:
			next(f)
			nonratios(w2, f, j, ratiolen+3)
			j += 1


workbook.close()
workbook2.close()

