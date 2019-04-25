#!/usr/bin/python3
# -*- coding: utf-8 -*-
#################################################################################################################
# Authors: Aleix Marine i Tena (aleix.marine@estudiants.urv.cat)                                                #
# Last review: 6/7/2018                                                                                         #
# Version: 1.0                                                                                                  #
#       																									    #
# Description: This code tries to predict interaction between genes in Pichia pastoris. We use supplementary    #
# materials from the article "Systematic analysis of complex genetic interactions"                              #
# (http://science.sciencemag.org/content/360/6386/eaao1729). DOI: 10.1126/science.aao1729 to get our data.      #
# We treat every gene as a node in our network. Links are every assay from the dataset between three genes.     #
# This links are tagged with 0 or 1 if there is interaction or not.                                             #
#                                                                                                               #
# Using Mixed-Membership Stochastic Block Model (MMSBM) we try to predict interaction between triplets of      	#
# genes. Every gene has a vector with the length of the number of groups. Every cell in this vector has the     #
# possibility of a gene behaving like one concrete group.                                                       #
#                                                                                                               #
# Permissions: Needs permission to read from data files and to write in the same folder as the script if        #
# to_File() method is called.                                                                                   #
#################################################################################################################
import getopt
import re
import codecs
import os
import random
import sys
import copy
import math
#import matplotlib.pyplot as plt
import numpy as np


# Contains all methods and data structures of the algorithm.
class Model:

	# Constructor:
	#
	# Description: Initializes data structures of algorithm.

	def __init__(self):
		# Vectors of probability, that relates the probability of a gene belonging to a determinate group
		self.ntheta = []
		self.theta = []

		# BIDIRECTIONAL DICTIONARY FOR USERS
		self.id_gene = {}  # dictionary that relates an id with its gene
		self.gene_id = {}  # dictionary that relates a gene with its id

		# Probability matrix
		self.pr = []
		self.npr = []

		# Matrix of ratings. rows: relation between gene1, gene2 and gene3 with the format "id1_id2_id3" in links,
		# and "namegene1_namegene2_namegene3" in nlinks.
		# columns: ratings (in this case 0 or 1). content: number of times seen relation between gene1, gene2
		#  and gene3 with rating r
		self.nlinks = {}
		self.links = {}

		# Same format as the previous dictionary, but we will use this data structure to calculate metrics, and not for
		# iterations
		self.test_links = {}

		# Contains the real and predicted probabilities for every triplete in the test_links
		self.results = []

		# Relates the id gene with its number of interactions
		self.uniqueg = {}

		# Describes how good our model explains our data
		self.likelihood = 0

		# Vector of likelihoods
		self.vlikelihood = []

		# R: number of possible ratings (tags) in a link. Typically 5 in a recommender. In our case, 1 for
		# interaction 0 for no interaction
		self.R = 2

		# K: Number of groups
		self.K = 0

		# P: Number of genes
		self.P = 0

		# eps: constant small value
		self.eps = 1e-10

	# Method initializeParameters:
	#
	# Description: Initializes theta and pr with random values and ntheta and npr with random values.
	# We will use pr and theta to do the iteration and npr and ntheta to store the new values
	#
	# Arguments:
	# 1.- K (number of group of genes) can be given using this method directly. Calling the method
	# without arguments will make K a random positive number
	#
	# Return parameters:
	# ntheta / theta 	--> vector of possibilities of a gene belonging to a determinate group of genes
	# npr / pr 			--> matrix of possibilities relating a triplet of group of genes in a certain rating
	# P 				--> num of genes
	# K 				--> number of groups of gene
	# R 				--> Number of possible ratings [0|1]

	def initialize_parameters(self, k_value=10):
		try:
			self.K = int(k_value)  # set K value to initialize parameters
		except ValueError:  # if we cant parse
			self.K = 10

		self.vlikelihood = []  # empty likelihood vector

		# Generate theta-like vectors
		self.theta = []
		self.ntheta = []
		for _ in range(self.P):  # Iterate over the number of genes
			a = [random.random() for _ in range(self.K)]  # xrange to generate big lists (optimization)
			self.theta.append(a)  # appends to theta a vector of random values
			self.ntheta.append([0.0] * self.K)  # generate a vector of reals with number of genes size and append it to ntheta

		# Generate pr and npr, 3D matrix with vectors of R components on its cells
		self.pr = []
		self.npr = []
		for i in range(self.K):
			b = []
			c = []
			for j in range(self.K):
				d = []
				e = []
				for k in range(self.K):
					a = [random.random() for _ in range(self.R)]
					d.append(a)
					e.append([0.] * self.R)
				b.append(d)
				c.append(e)
			self.pr.append(b)  # random values for pr
			self.npr.append(c)  # 0s for npr

		# Normalization for theta vector of genes:
		for i in range(self.P):  # Iterate over number of genes
			sumtheta = 0.  # sum of possibilities of theta vector for gene i
			for k in range(self.K):  # iterate over number of groups of genes
				sumtheta += self.theta[i][k]  # and get the sum of prob. of vector theta for gene 1

			if sumtheta < self.eps:  # Regenerate
					a = [random.random() for _ in range(self.K)]  # xrange to generate big lists (optimization)
					self.theta[i] = a

			sumtheta = sum(self.theta[i])
			for k in range(self.K):  # iterate over number of groups of genes
				try:
					self.theta[i][k] /= sumtheta
				except ZeroDivisionError:
					self.theta[i][k] /= (sumtheta + self.eps)  # adding an small value to avoid dividing by zero

		# Normalization of the vector probability for each gene having interaction with two other genes
		for i in range(self.K):
			for j in range(self.K):
				for k in range(self.K):
					sumpr = 0.  # Acumulator of possibilities for all ratings
					for r in range(self.R):
						sumpr += self.pr[i][j][k][r]
					for r in range(self.R):  # sum of prob of a user from group K giving rate R to a item from group L is 1
						try:
							self.pr[i][j][k][r] /= sumpr
						except ZeroDivisionError:
							self.pr[i][j][k][r] /= (sumpr + self.eps)

	# Method getInput:
	#
	# Description: Reads data from file in the same folder as this script, selects type of interaction and sets the
	# rating of the interaction.
	# All genes will be selected with P-value < 0.05.
	# P is also initialized with the number of genes.
	#
	# Arguments:
	# 1.- trainfile: File to read will be selected with this parameters
	# 2.- selectedInteractionType: [trigenic|digenic|*] Type of interaction to select. "*" to select all.
	# 3.- cutoffValue: selects the interaction as positive [1] with an adjusted genetic interaction score under
	# this value or as negative [0] if not. Use sys.float_info.max or sys.float_info.min in this parameter for assuming
	# all positive or negative interactions.
	# · digenic interaction cut-off in the original article (p < 0.05, |e| > 0.08)
	# · trigenic NOVEL interactions cut-off in the original article (p < 0.05, t < -0.08)
		# all trigenic novel and modified interactions p < 0.05 + trigenic
	# 4.- discard: [0|1] Add assays that are under the cutoffValue or discard them. By-default 0.
	# 5.- numlines: [0...sys.maxint] Allows to set the number of lines that are read from the dataset
	# Return parameters:
	# Fills with data the next variables: links, nlinks, id_gene, gene_id, P.
	#
	# File Format:
	# RAW Dataset S1
	# 1.- Query strain ID
	# 2.- Query allele name
	# 3.- Array strain ID
	# 4.- Array allele name
	# 5.- Combined mutant type
	# 6.- Raw genetic interaction score (epsilon)
	# 7.- Adjusted genetic interaction score (epsilon or tau)
	# 8.- P-value
	# 9.- Query single/double mutant fitness
	# 10.- Array single mutant fitness
	# 11.- Combined mutant fitness
	# 12.- Combined mutant fitness standard deviationç
	#
	# Treated Dataset S2
	# 1.- Query strain ID
	# 2.- Query allele name
	# 3.- Array strain ID
	# 4.- Array allele name
	# 5.- Combined mutant type
	# 6.- Adjusted genetic interaction
	# 7.- P-value
	# 8.- Interaction type

	def get_input(self, argfilename, selectedinteractiontype="trigenic", cutoffvalue=-0.08, discard=0, interactions='ALL'):
		try:
			gid = 0

			if selectedinteractiontype != 'trigenic' and selectedinteractiontype != 'digenic' and selectedinteractiontype != '*':
				raise ValueError("argument 2 selectedInteractionType must be trigenic, digenic or *")

			with codecs.open(argfilename, encoding='utf-8', mode='r') as fileref:

				line = fileref.readline()
				fields = re.split(r'\t+', line)  # Split line in fields using the regular expression '\t'
				if len(fields) == 12:
					reading_raw = 1
				else:
					reading_raw = 0

				for line in fileref.readlines():
					####
					# SELECT *
					# FROM dataset
					# WHERE P-value < 0.05
					# AND Combined mutant type = ["trigenic"|"digenic"|*]

					fields = re.split(r'\t+', line)  # obtain all fields from current line separated with tabs

					# dataset selection (now we can read from both types of dataset, s1 and s2)
					if reading_raw:
						fields.pop(5)

					# if current interaction type is not the type that we want to select, next element
					if selectedinteractiontype != "*":  # will be activated if selectedInteractionType different from *
						if fields[4] != selectedinteractiontype:
							continue

					# Decide if we want to consider positive vs negative interactions or novel vs digenic altered interactions
					if interactions == 'ALL':
						r = 0
						if float(fields[6]) < 0.05 and float(fields[5]) < cutoffvalue:
							r = 1  # check P-Value if < 0.05, R=1 - altered interactions
					else:
						# decide if positive or negative interaction taking in account cutoffValue
						if float(fields[6]) >= 0.05:  # check P-Value < 0.05, else, next
							continue

						if float(fields[5]) < cutoffvalue:
							r = 1
						else:
							if discard:
								continue
							else:
								r = 0

					# create list with the three implicated alleles
					gene_triplet = fields[1].split('+')
					gene_triplet.append(fields[3])

					# REGISTER ALLELES
					id_gene_triplet = []
					for gene in gene_triplet:  # for every gene in the triplet
						if gene not in self.gene_id.keys():  # gene hasn't been seen by algorithm
							self.gene_id[gene], n1 = gid, gid  # assign a new gid to this gene
							self.id_gene[gid] = gene  # assign a new gene to this gid
							self.uniqueg[n1] = 0  # when user identified by n1 is the first time found
							gid += 1  # update index gid
						else:  # gene is already registered
							n1 = self.gene_id[gene]  # get gid from gene, already registered

						# REGISTER NUMBER OF INTERACTIONS
						self.uniqueg[n1] += 1  # indicates number of interactions for gene identified by id

						# save ID from gene
						id_gene_triplet.append(str(n1))

					# Sort identifiers for unique key
					gene_triplet.sort()
					id_gene_triplet.sort()

					# Concatenate id and genes to create the key string
					str_name_gene_triplet = '_'.join(gene_triplet)
					str_gene_triplet = '_'.join(id_gene_triplet)  # joins genes with an underscore in between a triplet of genes

					try:
						self.links[str_gene_triplet][r] += 1  # link between g1, g2 and g3 with rating r it's been seen +1 times
						self.nlinks[str_name_gene_triplet][r] += 1
					except KeyError:  # if link between n1 and n2 with rating r is the first time seen then
						self.nlinks[str_name_gene_triplet] = [0] * 2
						self.nlinks[str_name_gene_triplet][r] += 1
						self.links[str_gene_triplet] = [0] * 2
						self.links[str_gene_triplet][r] += 1

				self.P = len(self.id_gene)  # get number of users
				fileref.close()

		except ValueError as error:
			print(error)
		except IOError as error:
			print('Error, file does not exist or can\'t be read')
			print(error)

	# Method that reads input links,nlinks and test_link form train and test files as 2name1_name2_name3\trating
	def get_traintest(self, trainfile, testfile):
		try:
			gid = 0

			with codecs.open(trainfile, encoding='utf-8', mode='r') as fileref:

				for line in fileref.readlines():

					fields = line.strip().split('\t')  # obtain all fields from current line separated with tabs

					gene_triplet = fields[0].split('_')
					r = int(fields[1])

					# REGISTER ALLELES
					id_gene_triplet = []
					for gene in gene_triplet:  # for every gene in the triplet
						if gene not in self.gene_id.keys():  # gene hasnt been seen by algorithm
							self.gene_id[gene], n1 = gid, gid  # assign a new gid to this gene
							self.id_gene[gid] = gene  # assign a new gene to this gid
							self.uniqueg[n1] = 0  # when user identified by n1 is the first time found
							gid += 1  # update index gid
						else:  # gene is already registered
							n1 = self.gene_id[gene]  # get gid from gene, already registered

						# REGISTER NUMBER OF INTERACTIONS
						self.uniqueg[n1] += 1  # indicates number of interactions for gene identified by id

						# save ID from gene
						id_gene_triplet.append(str(n1))

					# Sort identifier for unique key
					gene_triplet.sort()
					id_gene_triplet.sort()

					# Concatenate id and genes to create the key string
					str_name_gene_triplet = '_'.join(gene_triplet)
					str_gene_triplet = '_'.join(id_gene_triplet)  # joins genes with an underscore in between a triplet of genes

					try:
						self.links[str_gene_triplet][r] += 1  # link between g1, g2 and g3 with rating r it's been seen +1 times
						self.nlinks[str_name_gene_triplet][r] += 1
					except KeyError:  # if link between n1 and n2 with rating r is the first time seen then
						self.nlinks[str_name_gene_triplet] = [0] * 2
						self.nlinks[str_name_gene_triplet][r] += 1
						self.links[str_gene_triplet] = [0] * 2
						self.links[str_gene_triplet][r] += 1

				self.P = len(self.id_gene)  # get number of genes

			with codecs.open(testfile, encoding='utf-8', mode='r') as fileref:

				for line in fileref.readlines():

					fields = re.split(r'\t+', line)  # obtain all fields from current line separated with tabs

					gene_triplet = fields[0].split('_')
					r = int(fields[1])

					# REGISTER ALLELES
					id_gene_triplet = []
					for gene in gene_triplet:  # for every gene in the triplet
						if gene not in self.gene_id.keys():  # gene hasn't been seen by algorithm
							self.gene_id[gene], n1 = gid, gid  # assign a new gid to this gene
							self.id_gene[gid] = gene  # assign a new gene to this gid
							self.uniqueg[n1] = 0  # when user identified by n1 is the first time found
							gid += 1  # update index gid
						else:  # gene is already registered
							n1 = self.gene_id[gene]  # get gid from gene, already registered

						# REGISTER NUMBER OF INTERACTIONS
						self.uniqueg[n1] += 1  # indicates number of interactions for gene identified by id

						# save ID from gene
						id_gene_triplet.append(str(n1))

					# Sort identifiers for unique key
					gene_triplet.sort()
					id_gene_triplet.sort()

					# Concatenate id and genes to create the key string
					str_gene_triplet = '_'.join(id_gene_triplet)  # joins genes with an underscore in between a triplet of genes

					try:
						self.test_links[str_gene_triplet][r] += 1  # link between g1, g2 and g3 with rating r it's been seen +1 times
					except KeyError:  # if link between n1 and n2 with rating r is the first time seen then
						self.test_links[str_gene_triplet] = [0] * 2
						self.test_links[str_gene_triplet][r] += 1

				self.P = len(self.id_gene)  # get number of users

		except ValueError as error:
			print(error)
		except IOError as error:
			print('Error, file does not exist or can\'t be read')
			print(error)

		print('READ DATA train', len(self.links), len(self.nlinks))
		print('READ DATA test', len(self.test_links))

	# Method fold:
	#
	# Description: Produces test_sets and train_sets from the data_set contained in the model instance. Data that will
	# be split is contained in self.links, self.n_links and self.uniqueg.
	#
	# We will produce num_folds files of test set and num_folds files of train set. Every test set contains
	# (fraction * 100) % of the total triplets. There are no repeated triplets between test sets.
	# Every train set contains ( (1 - fraction) * 100) % of the total triplets. For every test set there is a corresponding
	# train test that fulfills the next equation train_test + test_set = whole_data_set (contained in self).
	#
	# There is a correspondence of a test set to a train set. We use a certain train_set to train the algorithm and the
	# test set to validate that training. All triplets of a certain test_set are not present in the corresponding train_
	# set.
	#
	# PRE: Data that is going to be split is stored in self.links and self.nlinks. self.uniqueg stores apparition of every
	# gene so is coherent with self.links and self.n_links dictionaries.
	# POST: All internal data structures remain UNTOUCHED. 2 * num_folds files will be created, each half for train or
	# test.
	#
	# Input Parameters:
	# - fraction: Relative size of test set

	def fold(self, fraction=0.2):
		test_set_size = int(len(self.links) * fraction)
		num_folds = int(1/fraction)

		# Add all triplets from self.links to array_links
		array_links = []
		for triplet, rating in self.links.items():
			array_links.append(triplet)
		np.random.shuffle(array_links)  # We shuffle the array_links list (every time we fold, output will be different)

		# Generate file pointers for all train and test files
		test_files = []
		train_files = []
		for i in range(num_folds):  # Create all train test file names & open file handles
			test_files.append(codecs.open('test' + str(i) + '.dat', encoding='utf-8', mode="w+"))
			train_files.append(codecs.open('train' + str(i) + '.dat', encoding='utf-8', mode="w+"))

		tests = []  # test stores num_folds splits of triplets
		# Append to test[] every slice of size test_set_size from array_links,
		for i in range(0, num_folds):
			tests.append(array_links[test_set_size * i: test_set_size * (i + 1)])

		# All remaining triplets in array links are put in the last test set
		tests[num_folds - 1] += array_links[test_set_size * num_folds:]

		# Output test files
		for num_fold in range(num_folds):
			for triplet in tests[num_fold]:
				rating = 1
				if self.links[triplet][0]:
					rating = 0
				ids = triplet.split("_")  # split triplet

				# Obtain names for the genes in the triplet.
				names = []
				for identifier in ids:
					name = self.id_gene[int(identifier)]  # obtain name
					names.append(name)  # accumulate

				names.sort()  # create identifier name string
				str_names = '_'.join(names)

				# file-print current test
				data = str_names + '\t' + str(rating) + '\n'
				test_files[num_fold].write(data)
				test_files[num_fold].flush()
			test_files[num_fold].close()

		# train stores num_folds training sets. Every training set contains all test_sets except the one that we are
		# testing against. train is a List of List of Strings.
		train = []
		for num_fold in range(num_folds):
			train_tmp = [] 	# Contains all triplets for a single train. train_tmp is a List of Strings
			for current_test in tests[: num_fold] + (tests[num_fold + 1:]):  # Select all tests except the counterpart
				train_tmp.extend(current_test)  # Acummulate all triplets from different tests in train_tmp
			train.append(train_tmp)
			for triplet in train[num_fold]:
				ids = triplet.split("_")  # split it
				# Obtain names for the genes in the triplet.
				names = []
				for identifier in ids:
					name = self.id_gene[int(identifier)]
					names.append(name)

				# obtain rating for the current triplet. (assuming just one given rating between r=0 or r=1)
				rating = 1
				if self.links[triplet][0]:
					rating = 0

				names.sort()
				str_names = '_'.join(names)

				# file-print current train
				data = str_names + '\t' + str(rating) + '\n'
				train_files[num_fold].write(data)
				train_files[num_fold].flush()
			train_files[num_fold].close()

	# Method do_prediction:
	#
	# Description: Returns the probability of interaction between three genes identified by IDs or names.
	#
	# Prerequisite: initialize_parameters, get_input, N x make_iteration. Needed to do predicitions.
	def do_prediction(self, id1, id2, id3):
		def calculate(i1, i2, i3):
			p = 0

			for i in range(self.K):
				for j in range(self.K):
					for k in range(self.K):
						p += self.theta[i1][i] * self.theta[i2][j] * self.theta[i3][k] * self.pr[i][j][k][1]  # We want positives ([1])
			return p

		try:
			id1_int, id2_int, id3_int = int(id1), int(id2), int(id3)
			probability = calculate(id1_int, id2_int, id3_int)
		except ValueError:
			probability = calculate(self.gene_id[id1], self.gene_id[id2], self.gene_id[id3])

		return probability

	# Method get_results:
	#
	# Description: Creates and returns the table of results from test set.
	#
	# Return Parameters: Results are an array of tuples of three elements, containing:
	# - The three genes involved in the interaction, coming from the tripletes in the test_set (self.test_links)
	# - Predicted probability of our model of this three genes interacting
	# - Real interaction [0|1]
	def calculate_test_set_results(self):
		self.results = []
		for triplet, rating in self.test_links.items():
			if rating[0]:
				rating = 0
			else:
				rating = 1
			id1, id2, id3 = triplet.split("_")
			self.results.append([self.do_prediction(id1, id2, id3), triplet, rating])

		# Sort results from less to higher predicted probability
		self.results.sort()
		self.results.reverse()

	# Method calculate_metrics:
	#
	# Description: Gets the cut_value taken as the threshold for choosing a sample as positive or negative.
	# I used a "rating" approach in which you compute the fraction of positives that are in your train_set.
	# Then you multiply this fraction by the number of samples in the test_set. You will obtain the number of predicted
	# positive in the train_set, assuming that both sets (train_set and test_set) are homogeneous distributed.
	#
	# We are going to sort our train_set by the prediction value, keeping higher prediction values on the top.
	# We will keep the first "number of predicted positives" values as positives, and the others as negatives.
	# We will use the predicted probability of interaction from the last element predicted as positiven as the cut_value.
	# This value will decide if a sample is positive (>= cut_value) or negative (< cut_value).

	def calculate_metrics(self):
		# calculate cut value using "ranking" method
		counter = 0
		cut_value = 0
		for triplet, rating in self.links.items():
			if rating[1] == 1:
				counter += 1
		positives_fraction = counter / len(self.links)  # Obtain the fraction of positives in our training set
		# Obtain the number of positives in the training set assuming that the distribution of 1 and 0 is equal between sets.
		positives_number = int(positives_fraction * len(self.test_links))
		counter = 0

		for data in self.results:
			if positives_number == counter:
				cut_value = data[0]
				break
			counter += 1

		# Calculate AUC metric
		positives = []
		negatives = []
		counter = 0
		for data in self.results:
			if data[2]:
				positives.append(data)
			else:
				negatives.append(data)

		for positive in positives:
			for negative in negatives:
				if positive[0] > negative[0]:
					counter += 1
		auc = counter / (len(positives) * len(negatives))

		# calculate metrics
		true_positives, false_positives, false_negatives, true_negatives = 0, 0, 0, 0
		for data in self.results:
			predicted = data[0]
			real = data[2]
			if predicted >= cut_value:
				if real:
					true_positives += 1
				else:
					false_positives += 1
			else:
				if real:
					false_negatives += 1
				else:
					true_negatives += 1

		precision = true_positives / (true_positives + false_positives)
		recall = true_positives / (true_positives + false_negatives)
		fallout = false_positives / (false_positives + true_negatives)

		return [precision, recall, fallout, auc]

	# Method get_data:
	#
	# Description: Reads data from a file indicated by argument with the to_string format, and copies it to
	# the different data structures in the model.
	def get_data(self, file_name="out0.txt"):

		def read_likelihood(file_reference):
			# Initialize likelihood vector
			self.likelihood = []
			while True:
				sam = str(file_reference.readline())
				if sam != '\n':
					sample_str, iterat_str, likelihood_str = sam.split('\t')
					try:
						samp, iterat, likelihood = int(sample_str), int(iterat_str), float(likelihood_str)
						self.likelihood.append([samp, iterat, likelihood])
					# If we can't parse, skip the line
					except ValueError:
						pass
				else:
					break
			return "read_likelihood"

		def read_list_of_genes(file_reference):
			# Initialize dictionaries
			self.gene_id = {}
			self.id_gene = {}
			self.uniqueg = {}
			while True:
				sam = str(file_reference.readline())
				if sam != '\n':
					gene_id_str, gene_name_str, numaparitions_str = sam.split('\t')
					try:
						gene_id, gene_name, numaparitions = int(gene_id_str), str(gene_name_str), int(numaparitions_str)
						self.gene_id[gene_name] = gene_id
						self.id_gene[gene_id] = gene_name
						self.uniqueg[gene_id] = numaparitions
					# If we can't parse, skip the line
					except ValueError:
						pass
				else:
					break
			return "read_list_of_genes"

		def read_links(file_reference):
			# Initialize dictionaries
			self.links = {}
			self.nlinks = {}
			while True:

				sam = str(file_reference.readline())
				if sam != '\n':
					# HARDCODED: Assuming R = 2
					gene_ids_str, aparitions1_str, aparitions0_str = sam.split('\t')
					try:
						gene_ids, aparitions1, aparitions0 = str(gene_ids_str), int(aparitions1_str), int(aparitions0_str)
						# Translate step
						id1, id2, id3 = gene_ids.split('\t')
						list_names = self.id_gene[id1], self.id_gene[id2], self.id_gene[id3],
						gene_names = '_'.join(list_names)

						try:
							self.nlinks[gene_names][0] = aparitions0
							self.nlinks[gene_names][1] = aparitions1
							self.links[gene_ids][0] = aparitions0
							self.links[gene_ids][1] = aparitions1
						except KeyError:
							self.links[gene_ids] = [0] * 2
							self.links[gene_ids][0] = aparitions0
							self.links[gene_ids][1] = aparitions1
							self.nlinks[gene_names] = [0] * 2
							self.nlinks[gene_names][0] = aparitions0
							self.nlinks[gene_names][1] = aparitions1
					# If we can't parse, skip the line
					except ValueError:
						pass
				else:
					break
			return "read_list_of_genes"

		def read_matrix(file_reference, matrix):
			for i in range(self.K):
				for _ in range(3):
					file_reference.readline()  # Skip three lines
				for j in range(self.K):
					sam = file_reference.readline()  # Get line
					list_data = sam.split('\t')  # Split line
					j = list_data.pop(0)
					for k in range(self.K):
						for r in range(self.R):
							matrix[i][j][k][r] = list_data.pop(0)
				for _ in range(2):
					file_reference.readline()  # skip two lines at the end
			return "read_list_of_genes"

		def read_vector(file_reference, vector_type):
			vector = []
			for _ in range(2):
				file_reference.readline()  # Skip the two first lines
			while True:
				data = file_reference.readline()
				if data != '\n':
					list_data = data.split('\t')
					id_gene = vector[list_data.pop(0)]
					vector[id_gene] = [] * self.K
					for k in range(self.K):
						vector[id_gene][k] = list_data.pop(0)
				else:
					break
			if vector_type == "theta":
				self.theta = vector
			else:
				self.ntheta = vector

		def read_rating(file_reference):
			self.R = int(file_reference.readline())

		def read_number_of_groups(file_reference):
			self.K = int(file_reference.readline())

		def switcher_method(argument):
			switcher = {
				"Number of groups of genes (K):": read_number_of_groups,
				"Number of possible ratings (R):": read_rating,
				"Likelihood vector:": read_likelihood,
				"LIST OF REGISTERED GENES": read_list_of_genes,
				"LIST OF LINKS BETWEEN GENE IDS": read_links,
				"MATRIX OF PROBABILITIES PR": lambda x: read_matrix(x, self.pr),
				"MATRIX OF PROBABILITIES NPR": lambda x: read_matrix(x, self.pr),
				"THETA VECTOR": lambda x: read_vector(x, self.theta),
				"NTHETA VECTOR": lambda x: read_vector(x, self.ntheta),
			}
			# Get the function from switcher dictionary
			return switcher.get(argument, lambda: argument+"is invalid line")

		self.initialize_parameters()
		try:
			file_ref = codecs.open(file_name, encoding='utf-8', mode="w+")
			for line in file_ref.readlines():
				lambda_function = switcher_method(line)
				# Execute the function
				lambda_function(file_ref)

			file_ref.close()

		except IOError:
			print("I/O error")

	# Method tostring:
	#
	# Description: Returns a CSV-like (using one tab as separator between fields) format string with data from the model
	# object.

	def to_string(self):

		def print_tuples(tuples):
			txt = ''
			txt += '\nPredicted Interaction\tID of genes\tReal Interaction\n'
			for tupla in tuples:
				txt += str(tupla[0])+'\t'+str(tupla[1])+'\t'+str(tupla[2])+'\n'
			return txt

		text = "Max Likelihood:\t" + str(self.likelihood) + "\n"
		text += "Number of genes (P):\t" + str(self.P) + "\n"
		text += "Number of links:\t" + str(len(self.links)) + "\n"
		text += "Number of groups of genes (K):\n" + str(self.K) + "\n"
		text += "Number of possible ratings (R):\n" + str(self.R) + "\n\n"
		self.calculate_test_set_results()  # Implicit call to calculate results
		text += "\nTest set:" + str(print_tuples(self.results))
		metrics = self.calculate_metrics()  # Implicit call to calculate metrics
		text += "\nMetrics:\nPrecision\tRecall\tFallout\tAUC\n"
		text += str(metrics[0]) + "\t" + str(metrics[1]) + "\t" + str(metrics[2]) + "\t" + str(metrics[3])

		# String of list of genes
		text += "\nLIST OF REGISTERED GENES\n"
		text += "Gene_ID\tGene_name\tnumAparitions\n"
		for gid in self.id_gene:
			text += str(gid) + "\t" + self.id_gene[gid] + "\t" + str(self.uniqueg[gid]) + '\n'

		return text

	# Method toFile(string):
	#
	# Description: Calls the toString function and prints it in a file. Overwrites the file if the file name given is
	# identical to the file name of a file in the same folder as this script.
	#
	# Arguments:
	# 1.- Name of the output file. By-default file name will be out.txt.

	def to_file(self, name_file=None):
		try:
			if name_file is None:
				name_file = "out.txt"
			fileref = codecs.open(name_file, encoding='utf-8', mode="w+")
			data = self.to_string()
			fileref.write(data)
			fileref.close()

		except IOError:
			print("I/O error")

	# Method compareLinks:
	#
	# Description: Given two object of type "Model" a and b, with a call "a.compareLinks(b)" we return a list of
	# links present in a that lack in b.
	#
	# Return parameters:
	# 1.- Returns "None" List if graph defined by "self" argument is a subgraph of "model" argument, otherwise, will
	# return the lacking links in b.

	def compare_links(self, arg_model):
		diff = []

		for link in self.nlinks.keys():
			if link in arg_model.nlinks.keys():
				continue
			else:
				diff.append(link)
		return diff

	# Method compareGenes:
	#
	# Description: Given two object of type "Model" a and b, with a call "a.compareLinksGenes(b)" we return a
	# list of genes present in a that lack in b.
	#
	# Return Parameters:
	# 1.- Returns "None" List if graph defined by "self" argument is a subgraph of "model" argument, otherwise, will
	# return the lacking genes in b.

	def compare_genes(self, arg_model):
		diff = []

		for gene in self.gene_id.keys():
			if gene in arg_model.gene_id.keys():
				continue
			else:
				diff.append(gene)
		return diff

	# Method computeLikelihood:
	#
	# Description: Computes likelihood of one of the current data sets contained in the class. By default the likelihood
	# of train set
	#
	# Return Parameters:
	# 1.- Returns the likelihood.
	def compute_likelihood(self, selectedset='train'):
		if selectedset == 'train':
			x = self.links.items()
		else:
			x = self.test_links.items()
		log_l = 0.
		for triplet, rating_vector in x:  # e X ra iteration
			g1, g2, g3 = triplet.split('_')
			id1, id2, id3 = int(g1), int(g2), int(g3)  # get IDs from three genes from link
			d = [self.eps] * self.R  # Generate a vector of R position with eps value (constant defined in headers)

			for i in range(self.K):
				for j in range(self.K):
					for k in range(self.K):
						for r in range(self.R):
							d[r] += self.theta[id1][i] * self.theta[id2][j] * self.theta[id3][k] * self.pr[i][j][k][r]
			for r in range(self.R):
				log_l += rating_vector[r] * math.log(d[r])
		self.likelihood = log_l
		return log_l

	# Method makeIteration:
	#
	# Description: Do recursive computation to iterate over the pr, and theta vector. New values will be stored in
	# npr and ntheta instance variables. Then, new values are normalized.
	#
	# Return Parameters:
	# Updates and normalizes new values in npr and ntheta data structures.

	def make_iteration(self):

		counter = [0] * self.P
		for triplet, rating_vector in self.links.items():  # e X R iteration
			g1, g2, g3 = triplet.split('_')
			id1, id2, id3 = int(g1), int(g2), int(g3)  # get IDs from three genes from link
			d = [self.eps] * self.R  # Generate a vector of R position with eps value (constant defined in headers)

			counter[id1] += 1
			counter[id2] += 1
			counter[id3] += 1

			for i in range(self.K):
				for j in range(self.K):
					for k in range(self.K):
						for r in range(self.R):
							d[r] += self.theta[id1][i] * self.theta[id2][j] * self.theta[id3][k] * self.pr[i][j][k][r]

			for i in range(self.K):
				for j in range(self.K):
					for k in range(self.K):
						for r in range(self.R):
							# auxiliary variable
							a = (self.theta[id1][i] * self.theta[id2][j] * self.theta[id3][k] * self.pr[i][j][k][r]) / d[r]
							self.ntheta[id1][i] += a * rating_vector[r]
							self.ntheta[id2][j] += a * rating_vector[r]
							self.ntheta[id3][k] += a * rating_vector[r]
							self.npr[i][j][k][r] += a * rating_vector[r]

		# Normalizations:
		# divide all possibilities of i belonging to a group k with the number of relation of that user
		for i in range(self.P):
			for k in range(self.K):
				self.ntheta[i][k] /= float(counter[i])

		# divide the probability of the group k giving rate to a item l with a rating r between the sum of all ratings
		for i in range(self.K):
			for j in range(self.K):
				for k in range(self.K):
					d = self.eps
					for r in range(self.R):
						d += self.npr[i][j][k][r]
					for r in range(self.R):
						self.npr[i][j][k][r] /= d

		# Copies values from n* data structures to the current structures
		self.theta = copy.copy(self.ntheta)
		for i in range(self.K):
			for j in range(self.K):
				for k in range(self.K):
					self.pr[i][j][k] = self.npr[i][j][k]

		# Reinitialization of n * data structures
		for i in range(self.P):
			self.ntheta[i] = [0.] * self.K
		for i in range(self.K):
			for j in range(self.K):
				for k in range(self.K):
					self.npr[i][j][k] = [0.] * self.R

	# Method compareDataset(model):
	#
	# Description: Given another object model, data of links and genes is comparated separately. This information will
	# be printed in the output.
	#
	# Return parameters:
	# 1.- True if the first dataset is subgraph for genes and links, False if not.

	def compare_dataset(self, arg_model):
		if not self.compare_links(arg_model):
			print("First dataset is subgraph of second dataset for links")
			node = 1
		else:
			print("First dataset is not subgraph of second dataset for links")
			node = 0
		if not self.compare_genes(arg_model):
			print("First dataset is subgraph of second dataset for nodes")
			link = 1
		else:
			print("First dataset is not subgraph of second dataset for nodes")
			link = 0

		return link and node


# Function compareS1withS2:
#
# Description: Checks that the filtered dataset and the raw dataset are equal. The code in this function is hardcoded
# and its unique purpose is demonstrate that the selection criteria specified in article is consistent
# with obtained data from filtering dataset S1 for trigenic interaction.


def compares1withs2():
	rawmodel = Model()
	rawmodel.get_input('Data_S1.csv', 'trigenic', -0.08, 1)  # discards negatives for raw dataset

	treatedmodel = Model()
	treatedmodel.get_input('Data_S2.csv', 'trigenic', sys.float_info.max)  # take all data from treated dataset

	treatedmodel.to_file("treated.txt")
	rawmodel.to_file("raw.txt")

	print("\nComparing treated with raw: ")
	if treatedmodel.compare_dataset(rawmodel):
		sub0 = 1
	else:
		sub0 = 0

	print("\nComparing raw with treated: ")
	if rawmodel.compare_dataset(treatedmodel):
		sub1 = 1
	else:
		sub1 = 0

	if sub1 and sub0:
		print("\nDatasets are equal")
	else:
		print("\nDatasets are different")


# help Function:
#
# Description: Gives a small help about the usage of the algorithm

def usage(it, s, check, file, k):
	txt = '\n\nUsage:\n'
	txt += './TrigenicInteractionPredictor.py [-h|--help][-i|--iterations=] <number_of_iterations> '
	txt += '[-s|--samples=] <number_of_samples> [-c|--check=] <frequency check> [-f|--file=] <file name> '
	txt += '[-k|--k=] <Number of groups>'
	txt += '\n\n\nDescription of the arguments:\n[-h|--help] Calls help and exits the program.\n[-i|--iterations=]'
	txt += ' Number of iterations done per sample for training the algorithm. More iterations increase the '
	txt += 'likelihood of the model.\n[-s|--samples=] Number of times the algorithm is computated.\n[-c|--check=]'
	txt += ' Sets the number of iterations that need to pass for calculating the likelihood. Set it to 0 for just'
	txt += ' checking the likelihood after computating every sample. Calculate the likelihood frequently decreases'
	txt += ' the throughput of the algorithm.\n[-f|--file=] Relative (from this file) or absolute path to the'
	txt += ' Data File.\n[-t|--type=] Selects which kind of interactions are going to be selected. Three possible'
	txt += ' options available, trigenic, digenic or all, for selecting all types of interactions. In version 1.0, the'
	txt += ' algorithm is not capable of processing digenic interactions.\n[-v|--cutvalue=] Value used for '
	txt += 'determining if an interaction is positive or negative. \n[-k|--k=] Number of groups.\n\nArguments '
	txt += 'can be left blank and the next default values will be selected:\n Number of iterations per sample: ' + str(
		it)
	txt += '\nNumber of samples: ' + str(s) + '\nLikelihood frequency checking: ' + str(
		check) + '\nData input file: ' + str(file)
	txt += '\nNumber of groups: ' + str(k)
	txt += '\n\nThis code is optimized to be executed by pypy3. You can find a pypy3 environment in'
	txt += 'the src folder of the project. In a terminal situated in that folder, the execution with pypy3'
	txt += 'follows the next pattern:\n\n./pypy3/bin/pypy3 TrigenicInteractionPredictor.py {arguments}\n'
	print(txt)


# Main Function:
#
# Description: Computes the algorithm given a dataset and a number of iterations.
#
# Arguments:
# 1.- iterations: Number of iterations done by algorithm.
# 2.- samples: Number of samples done by algorithm.
# 3.- frequencyCheck: Number of iterations needed to check if likelihood has converged.
# 4.- trainfile: Name of the dataset trainfile.
# 5.- interactionType: Type of interaction selected.
# 6.- cutOffValue: Value used to determine if an interaction is positive or negative.


if __name__ == "__main__":
	random.seed(os.getpid())

	# Default arguments
	iterations = 10000
	samples = 100
	frequencyCheck = 10
	train = ""
	test = ""
	argk = 4

	# This method returns value consisting of two elements: the first is a list of (option, value) pairs.
	# The second is the list of program arguments left after the option list was stripped.
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:s:c:t:e:k:", ["help", "iterations=", "samples=", "check=", "train=", "test=", "k="])

		for opt, arg in opts:
			if opt in ("-h", "--help"):  # Show the usage if help is called
				usage(iterations, samples, frequencyCheck, train, argk)
				exit(0)
			elif opt in ("-i", "--iterations"):
				if int(arg) < 1:
					print("\n\nERROR: Number of iterations should be a integer positive number!")
					raise ValueError
				iterations = int(arg)
			elif opt in ("-s", "--samples"):
				if int(arg) < 1:
					print("\n\nERROR: Number of samples should be a integer positive number!")
					raise ValueError
				samples = int(arg)
			elif opt in ("-c", "--check"):
				if int(arg) < 0:
					print("\n\nERROR: frequency of checking should be a integer positive number or 0!")
					raise ValueError
				if int(arg) == 0:
					frequencyCheck = iterations + 1  # Checking likelihood just at the end
				else:
					frequencyCheck = int(arg)
			elif opt in ("-t", "--train"):
				if os.path.isfile(arg):
					train = arg
				else:
					print("\n\nERROR: The selected file does not exist.")
					raise ValueError
			elif opt in ("-e", "--test"):
				if os.path.isfile(arg):
					test = arg
				else:
					print("\n\nERROR: The selected file does not exist.")
					raise ValueError
			elif opt in ("-k", "--k"):
				if int(arg) < 1:
					print("\n\nERROR: Number of groups should be a positive integer number different from 0")
					raise ValueError
				if int(arg == 1):
					print(
						"\n\nWARNING:If number of groups is 1, algorithm is single-membership instead of mixed\n")
				argk = int(arg)

	except getopt.GetoptError:
		usage(iterations, samples, frequencyCheck, train, argk)
		sys.exit(2)
	except ValueError:
		usage(iterations, samples, frequencyCheck, train, argk)
		sys.exit(2)

	# Display Warning
	if int(frequencyCheck) > int(iterations):
		print("\n\nWARNING: the likelihood frequency checking is bigger that the number of iterations per sample.")
		print("likelihood will only be calculated at the end of every sample (equivalent to --check=0 or -c 0)\n")

	# Start Algorithm
	msg = "\n****************************************\n* Trigenic Interaction Predictor v 1.0 *\n**************"
	msg += "**************************\n\nDoing " + str(samples) + " samples of " + str(iterations) + " iterations."
	msg += "\nData is read from train-file " + train + " and test-file" + test + ".\n"
	msg += "K value (number of groups) is " + str(argk) + "."
	msg += "\nLikelihood will be calculated every " + str(frequencyCheck) + " iterations."
	print(msg)

	model = Model()
	model.get_traintest(train, test)

	print("\nStarting algorithm...")

	for sample in range(int(samples)):
		print("Sample " + str(1 + sample) + ":")
		model.initialize_parameters(argk)
		print("Parameters have been initialized")
		like0 = model.compute_likelihood()
		print("·Initial Likelihood is "+str(like0))
		model.vlikelihood.append([sample, 0, like0])  # append result into the global vector of likelihoods

		for iteration in range(iterations):
			model.make_iteration()

			if iteration % frequencyCheck == 0:
				like = model.compute_likelihood()
				print("· Likelihood " + str(iteration + 1) + " is " + str(like))
				model.vlikelihood.append([sample, iteration + 1, like])  # append result into the global vector of likelihoods
				if math.fabs((like - like0) / like0) < 0.0001:
					print("\n\t****************************\n\t* Likelihood has converged *\n\t****************************")
					outfile = 'Sample' + str(sample) + '_K' + str(argk) + '.csv'
					model.to_file(outfile)  # // debug

					break
				like0 = like
