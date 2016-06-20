import sys
import xlsxwriter
import os

#python table.py "scores+positions(1-2-2-1).txt" "scores+positions(2-2-5-4).txt" "scores+positions(3-2-4-4).txt" "scores+positions(5-1-8-4).txt" "scores+positions(5-2-8-3).txt" "scores+positions(5-3-8-2).txt" "scores+positions(5-4-6-4).txt"

def firstcol(w, f):
	row = 1
	col = 0
	
	for line in f:
		line = line.split(" ")		
		
		w.write(row, col, line[0]+" - "+line[1])
		
		row += 1 
	

def firstline(k, w, name):	

	w.write(0, k, "Alignment Score " +str(k)+": "+ name[16:-4] + " / Position")

			

def writecol(i, f, w):
	print i
	row = 1
	col = i
	
	scores = []
		
	for line in f:
		if line == "":
			break
				
		line = line[:-1].split(" ")
		
		scores.append(line)
		
	found = False 
	count = 0
	countall = 0
	
		
	for s in scores:
		
		countall = len(scores)
	
		posis = []
		endstr = ""
		 
		key = s[0] + " " + s[1]
		
		endstr += str(s[2])+" "+str(s[3])+" "
		
		if key in pos:
			posis = pos[key]	
			
			if len(posis) == 0:
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
		found = False	
		

		
#main

workbook = xlsxwriter.Workbook('func_result_table.xlsx')
w = workbook.add_worksheet()

w.write(0, 0, "MTI")

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



for i in range (0, len(filelist)):
	
	if first:		
		if filelist[i][0:6] == "scores":		
			with open(str(path+"\\"+filelist[i])) as one:
				firstcol(w, one)
			
			first = False
		else:
			continue
	
	if filelist[i][0:6] == "scores":	
		with open(path+"\\"+filelist[i]) as f:	
			firstline(k, w, str(filelist[i]))	
			writecol(k, f, w)
			k += 1
			
		

workbook.close()

