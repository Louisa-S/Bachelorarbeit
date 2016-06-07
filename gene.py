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
		#print genelist[0]
		return genelist
