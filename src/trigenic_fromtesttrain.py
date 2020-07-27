#!/usr/bin/python3
# -*- coding: utf-8 -*-
#################################################################################################################
## Main that uses the code in TrigenicInteractionPredictor but reading from existing
## training and test files using the get_traintest method                                                       #
#################################################################################################################

import re
import codecs
import os
import random
import sys
import copy
import math
#import matplotlib.pyplot as plt
import numpy as np
from TrigenicInteractionPredictor import *
# Main Function:
#
# Description: Computes the algorithm given a dataset and a number of num_iterations.
#
# Arguments:
# 1.- num_iterations: Number of num_iterations done by algorithm.
# 2.- samples: Number of samples done by algorithm.
# 3.- frequency_check: Number of num_iterations needed to check if likelihood has converged.
# 4.- filename: Name of the dataset filename.
# 5.- interactionType: Type of interaction selected.
# 6.- cutOffValue: Value used to determine if an interaction is positive or negative.


if __name__ == "__main__":
	random.seed(os.getpid())
        

	# READ ARGUMENTS
	try:
		iterations = int(sys.argv[1])
		samples = int(sys.argv[2])
		frequencyCheck = int(sys.argv[3])
		argk = int(sys.argv[4])
		sampleini = int(sys.argv[5])
		trainfile = sys.argv[6]
		testfile = sys.argv[7]
                
                
	# BY-DEFAULT VALUES
	except IndexError:
		iterations = 10000
		samples = 100
		frequencyCheck = 10
		argk = 10
		itini = 0
		trainfile = 'train.dat'
		testfile = 'test.dat'
	
	msg = "\n****************************************\n* Trigenic Interaction Predictor v 1.0 *\n**************"
	msg += "**************************\n\nDoing "+str(samples)+" samples of "+str(iterations)+" num_iterations."
	msg += "**************************\n\nStarting from sample "+str(sampleini)+" ."
	msg += "\nLikelihood will be calculated every "+str(frequencyCheck)+" num_iterations."
	print(msg)

	model = Model()
	model.get_train_test(trainfile, testfile)
#	model.fast_fold(output=1)

	print("\nStarting algorithm...")

	for sample in range(sampleini,sampleini+int(samples)):
		print("Sample "+str(1 + sample)+":")
		model.initialize_parameters(argk)
		print("Parameters have been initialized")
		like0 = model.compute_likelihood()
		print("· Likelihood 0 is "+str(like0))
		model.vlikelihood.append([sample, 0, like0])  # append result into the global vector of likelihoods

		for iteration in range(iterations):
			model.make_iteration()

			if iteration % frequencyCheck == 0:
				like = model.compute_likelihood()
				print("· Likelihood " + str(iteration + 1) + " is " + str(like))
				model.vlikelihood.append([sample, iteration + 1, like])  # append result into the global vector of likelihoods
				if math.fabs((like - like0) / like0) < 0.001:
					print("\n\t****************************\n\t* Likelihood has converged *\n\t****************************")

#					model.to_file("outprev" + str(sample) + ".csv") # // debug
#					model.get_data("out1.txt")
					outfile = 'outSamp%dK%d.csv' % (sample,argk)
#					model.to_file("late"+str(sample)+".csv") # // debug
					model.to_file_short(outfile) # // debug
#					model.calculate_test_set_results()
#					precision,recall, fallout,auc = model.calculate_metrics()
#					print('Precision,',precision)
#					print('Recall,',recall)
#					print('AUC,',auc)
                                        
					break
				like0 = like
			# model.get_data("out2.txt")  # // debug
		#model.to_file("out" + str(sample) + ".csv")
