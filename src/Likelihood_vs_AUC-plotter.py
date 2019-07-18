#!/usr/bin/python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import codecs
import re
import numpy as np
import os

if __name__ == "__main__":

	# Create data
	#os.system(more /media/sf_DEFINITIVE_RESULTS/K*/fold*/Sample* | grep -oP "fold[0-5]/Sample_[0-9]+_K.\.csv|(Held-out Likelihood:\t-[0-9]*\.[0-9]*)|^([0-9]\.[0-9]*\t){3}[0-9]\.[0-9]*")
	results = []

	for k in range(4):
		results.append([[[], []], [[],[]], [[], []], [[], []], [[], []]])

	with codecs.open('/home/aleixmt/Desktop/TrigenicInteractionPredictor/data/Filtered_results.csv', encoding='utf-8', mode='r') as fileref:
		while True:
			line = fileref.readline()
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




'''
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
'''