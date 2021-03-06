#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import codecs

resultsFolder = "../data/REDUCED_TRAIN_TEST_RESULTS/"

results = []
for i in range(4):  # List of 4 positions corresponding to each K //HC
	a = []
	for j in range(5):  # List of 5 positions corresponding to each fold.
		b = []  # Will be a list of 2 positions, likelihood and AUC
		a.append(b)
	results.append(a)

for k in range(4):
	for fold in range(5):
		with codecs.open(resultsFolder + "K" + str(int(k) + 2) + "_fold" + str(fold) + ".csv", encoding='utf-8', mode='r') as fileHandler:
			fileHandler.readline()
			fileHandler.readline()
			line = fileHandler.readline()
			fields = line.split('\t')
			HeldOutLikelihood = float(fields[0])
			AUC = float(fields[1])
			results[int(k)][fold].append(AUC)
			results[int(k)][fold].append(HeldOutLikelihood)

# Create plot
fig, ax = plt.subplots()

colors = ('black', 'red', 'blue', 'purple')
markers = ('.', 's', 'p', 'P', 'd')

for k in range(4):
	x, y = results[k][0]
	ax.scatter(x, y, c=colors[k], marker=markers[0], label="K" + str(int(k) + 2), s=40)

for fold in range(5):
	x, y = results[0][fold]
	ax.scatter(x, y, c=colors[0], marker=markers[fold], label="fold " + str(fold), s=40)

for k in range(1, 4):
	for fold in range(1, 5):
		x, y = results[k][fold]
		if fold == 0:
			ax.scatter(x, y, c=colors[k], marker=markers[fold], s=40)
		elif k == 0:
			ax.scatter(x, y, c=colors[k], marker=markers[fold], s=40)
		else:
			ax.scatter(x, y, c=colors[k], marker=markers[fold], s=40)



plt.title('Likelihood vs AUC')
plt.legend()
plt.grid(True)
plt.xlabel("AUC")
plt.ylabel("Likelihood")

plt.show()
'''
				if re.match("Sample_[0-9]*", str(line)):
					pass
				else:
					break
				# Decide where the current result belongs
				K_index = -1
				if re.search("K5", line):
					K_index = 3
				elif re.search("K4", line):
					K_index = 2
				elif re.search("K3", line):
					K_index = 1
				elif re.search("K2", line):
					K_index = 0

				if re.search("fold0", line):
					fold_index = 0
				elif re.search("fold1", line):
					fold_index = 1
				elif re.search("fold2", line):
					fold_index = 2
				elif re.search("fold3", line):
					fold_index = 3
				elif re.search("fold4", line):
					fold_index = 4

				line = fileref.readline()  # Obtain line of the likelihood
				likelihood = line.split("\t")[1]

				metrics = fileref.readline()  # Obtain AUC metric
				AUC = metrics.split("\t")[3].split("\n")[0]

				results[K_index][fold_index][0].append(float(likelihood))
				results[K_index][fold_index][1].append(float(AUC))

print("empiesa la funsion")
# Create plot
fig, ax = plt.subplots()
i = 0
for color in ['red', 'green', 'blue', 'yellow']:
	y, x = data[i]
	scale = 10.0 * 0.1
	ax.scatter(x, y, c=color, label="K"+str(i+2), s=4)
	i = i + 1

plt.title('Likelihood vs AUC')
plt.legend()
plt.grid(True)
plt.xlabel("AUC")
plt.ylabel("Likelihood")

plt.show()




fold4 es el que dona mes aalt
Fold validatipon per a triar k
Kold 3 aparentment 

Overleaf: editoonline de latex ( amb plantillles)
texmaker

Agafar scores i fer miitjana / mediana

Plot de la mitjana de cada fold

Fer mitjana de cada probbabilitat a cada sample de cada tripleta testejada. 
Obtenir una metrica de mediana mitjana per a cada trippeta per a recalcular AUC i held-out likelihood (mitjana)
Representar cada held-out likelihood mig i AUC mig front K 


Gene ontology (gene ontology slim term mapper)

Biopython 

'''
