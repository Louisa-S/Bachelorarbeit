import urllib2
import sys

def mirid(f):
	ids = dict()
	
	for line in f:
	
		
		if "Non-Functional" in line:
			continue
		if "Functional MTI (Weak)" in line:
			continue
			
		line = line.split(";")
		
		
		
		if line[0] not in ids:
			ids[line[0]] = [line[1], line[4]]
	
	
	return ids
		

def parseurl(mlist):
	
	positions = dict()
	
	for m in mlist:	
		
		posis = []	

		url = 'http://mirtarbase.mbc.nctu.edu.tw/php/detail.php?mirtid='+m+'#target'
		
		data = urllib2.urlopen(url)
		
			
		f = data.readlines()

		s = []
		schreiben = False

		for line in f:
			#print line
			
			if "<td>NM_" in line:
				 l1 = line.replace(" ", "")
				 l1 = l1.split("&")
				 posis.append(l1[0][4:])				 
				 
			
			if "<th>miRNA-target interactions (Pred" in line:
				schreiben = True
	
			if schreiben:
				if "</table>" in line:
					break
				s.append(line)
			else:
				continue		
		
		for line in s:
			if " - " in line:
				line = line.replace(" ", "")
				pos = line.split("-")
				posis.append(pos[0])
		
		if posis == []:
			continue
		
		positions[mlist[m][0]] = posis
		
		print posis
		
		
				
	return positions
		
		
with open(str(sys.argv[1])) as f:
	next(f)
	mlist = mirid(f)
	
	#print mlist
		
	positions = parseurl(mlist)
	
	res = positions.items()
	
	with open("mirtargets_positions.txt", "w") as out:
		for i in res:
			out.write(str(i[0]+" "))
			
			for k in range(len(i[1])):
				out.write(str(i[1][k])+" ")
			
			out.write("\n")

























