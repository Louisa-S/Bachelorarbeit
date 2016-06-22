import sys
import matplotlib.pyplot as plt

tar = []
non = []

with open("./results/scores+positions5-3-8-2.txt") as target:
	for line in target:
		if line == "\n":
			break
		else:
			line = line.split()
			tar.append(int(line[2][:-2]))

with open("./results/non-scores+positions5-3-8-2.txt") as nontarget:
	for line in nontarget:
		if line == "\n":
			break
		else:
			line = line.split()
			non.append(int(line[2][:-2]))


tar = sorted(tar)
non = sorted(non)

scores = dict()
nonscores = dict()


for t in tar:
	if t in scores:
		scores[t] += 1
	else:
		scores[t] = 1

for n in non:
	if n in nonscores:
		nonscores[n] += 1
	else:
		nonscores[n] = 1

p = []
plotx = []
ploty = []

print scores
print nonscores

allscores = []

for key in scores:
	if key not in allscores:
		allscores.append(key)
		
for key in nonscores:
	if key not in allscores:
		allscores.append(key)
		

for key in allscores:
	if key in scores:
		x = float(scores[key])
	else:
		x = 0
	if key in nonscores:
		y = float(nonscores[key])
	else:
		y = 0
	
	ratio = float(x) / float(x+y)
	
	p.append((key, ratio))
	
	
for part in p:
	print part
	plotx.append(part[0])
	ploty.append(part[1])
	
	

plt.bar(plotx, ploty)
plt.xlabel("Alignment score")
plt.ylabel("Percentages of real targets")
plt.savefig("results/plot_scores.png")
