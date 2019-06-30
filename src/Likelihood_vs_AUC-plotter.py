#!/usr/bin/python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import codecs
import re
import numpy as np
import os

if __name__ == "__main__":

	# Create data

	results = []
	K2 = [[], []]
	K3 = [[], []]
	K4 = [[], []]
	K5 = [[], []]

	results.append(K2)
	results.append(K3)
	results.append(K4)
	results.append(K5)

	with codecs.open("/home/aleixmt/Escritorio/TrigenicInteractionPredictor/data/Filtered_results.csv", encoding='utf-8', mode='r') as fileref:
		while True:
			# Decide where the current result belongs
			line = fileref.readline()

			if re.match("K.*", line):  # shitty way to know if we are at the end (fileread returns garbage)
				pass
			else:
				break

			K_index = -1
			if re.search("K5", line):
				K_index = 3
			elif re.search("K4", line):
				K_index = 2
			elif re.search("K3", line):
				K_index = 1
			elif re.search("K2", line):
				K_index = 0

			line = fileref.readline()  # Obtain line of the likelihood
			likelihood = line.split("\t")[1].split("\n")[0]


			metrics = fileref.readline()  # Obtain AUC metric
			AUC = metrics.split("\t")[3].split("\n")[0]

			results[K_index][0].append(float(likelihood))
			results[K_index][1].append(float(AUC))

	print("empiesa la funsion")
	data = (K2, K3, K4, K5)
	# Create plot
	fig, ax = plt.subplots()
	i = 0
	for color in ['red', 'green', 'blue', 'yellow']:
		y, x = data[i]
		scale = 200.0 * 0.6
		ax.scatter(x, y, c=color, label="K"+str(i+2))
		i = i + 1

	plt.title('Likelihood vs AUC')
	plt.legend()
	plt.grid(True)
	plt.xlabel("AUC")
	plt.ylabel("Likelihood")

	plt.show()



