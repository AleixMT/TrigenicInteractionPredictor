#!/usr/bin/python
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
# Using Mixed-Membership Stochastic Block Model (MMSBM) we try to predict interaction between tripletes of      #
# genes. Every gene has a vector with the length of the number of groups. Every cell in this vector has the     #
# possibility of a gene behaving like one concrete group.                                                       #
#                                                                                                               #
# Permissions: Needs permission to read from data files and to write in the same folder as the script if        #
# toFile() method is called.                                                                                    #
#################################################################################################################

import re
import codecs
import os
import random
import sys
import time
import copy
import math


class Model:

	# Constructor:
	#
	# Description: Initializes data structures of algorithm.

	def __init__(self):
		# Vectors of probability, that relates the probability of user u belonging to a determinate group
		self.ntheta = []
		self.theta = []

		# BIDIRECTIONAL DICTIONARY FOR USERS
		self.id_gene = {}  # dictionary that relates an id with its gene
		self.gene_id = {}  # dictionary that relates a gene with its id

		# Probability matrix
		self.pr = []
		self.npr = []

		# Matrix of ratings. rows: relation between gene1, gene2 and gene3 wit the format "id1_id2_id3".
		# columns: ratings (in this case 0 or 1). content: number of times seen relation between gene1, gene2
		#  and gene3 with rating r
		self.links = {}

		# Matrix of ratings. rows: relation between gene1, gene2 and gene3 wit the format "gene1_gene2_gene3".
		# columns: ratings (in this case 0 or 1). content: number of times seen relation between gene1, gene2
		#  and gene3 with rating r
		self.nlinks = {}

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
	# Description: Initializes theta and pr with random values and ntheta and npr with random values
	#
	# Arguments:
	# 1.- K (number of group of genes) can be given using this method directly. Calling the method
	# without arguments will make K a random positive number
	#
	# Return parameters:
	# ntheta / theta 	--> vector of possibilities of a gene belonging to a determinate group of genes
	# npr / pr 			--> matrix of possibilities relating a triplete of group of genes in a certain rating
	# P 				--> num of genes
	# K 				--> number of groups of gene
	# R 				--> Number of possible ratings [0|1]

	def initializeParameters(self, k=None):

		if k is None:
			k = 0
		while k <= 0:  # if k is 0 or negative initialize k with a random value
			k = random.random()

		self.setK(int(k))  # set K value to initilize parameters

		for _ in range(self.P):  # Iterate over the number of genes
			a = [random.random() for _ in xrange(self.K)]  # xrange to generate big lists (optimization)
			self.theta.append(a)  # appends to theta a vector of random values
			self.ntheta.append([0.0] * self.K) # generate a vector of reals with number of genes size and append it to ntheta

		# Generate pr and npr, 3D matrix with vectors of R components on its cells
		for i in range(self.K):
			b = []
			c = []
			for j in range(self.K):
				d = []
				e = []
				for k in range(self.K):
					a = [random.random() for _ in xrange(self.R)]
					d.append(a)
					e.append([0.] * self.R)
				b.append(d)
				c.append(e)
			self.pr.append(b)  # random values for pr
			self.npr.append(c)  # 0s for npr

		# Normalization for theta vector of genes:
		for i in range(self.P):  # Iterate over number of genes
			sumTheta = 0.  # sum of possibilities of theta vector for gene i
			for k in range(self.K):  # iterate over number of groups of genes
				sumTheta += self.theta[i][k]  # and get the sum of prob. of vector theta for gene 1
			for k in range(self.K):  # normalization for all components,
				try:
					self.theta[i][k] /= sumTheta
				except ZeroDivisionError:
					self.theta[i][k] /= (sumTheta + self.eps)  # adding an small value to avoid dividing by zero

		# Normalization of the vector probability for each gene having interaction with two other genes
		for i in range(self.K):
			for j in range(self.K):
				for k in range(self.K):
					sumPr = 0.  # Acumulator of possibilities for all ratings
					for r in range(self.R):
						sumPr += self.pr[i][j][k][r]
					for r in range(self.R):  # sum of prob of a user from group K giving rate R to a item from group L is 1
						try:
							self.pr[i][j][k][r] /= sumPr
						except ZeroDivisionError:
							self.pr[i][j][k][r] /= (sumPr + self.eps)

	# Method getInput:
	#
	# Description: Reads data from file in the same folder as this script, selects type of interaction and sets the
	# rating of the interaction.
	# All genes will be selected with P-value < 0.05.
	# P is also initialized with the number of genes.
	#
	# Arguments:
	# 1.- filename: File to read will be selected with this parameters
	# 2.- selectedInteractionType: [trigenic|digenic|*] Type of interaction to select. "*" to select all.
	# 3.- cutoffValue: selects the interaction as positive [1] with an adjusted genetic interaction score under
	# this value or as negative [0] if not. Use sys.float_info.max or sys.float_info.min in this parameter for assuming
	# all positive or negative interactions.
	# · digenic interaction cut-off in the original article (p < 0.05, |e| > 0.08)
	# · trigenic interactions cut-off in the original article (p < 0.05, t < -0.08)
	# 4.- discard: [0|1] Add assays that are under the cutoffValue or discard them. By-default 0.
	#
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

	def getInput(self, filename, selectedInteractionType, cutoffValue, discard=0):
		try:
			gid = 0

			if selectedInteractionType != 'trigenic' and selectedInteractionType != 'digenic' and selectedInteractionType != '*':
				raise ValueError("argument 2 selectedInteractionType must be trigenic, digenic or *")

			with codecs.open(filename, encoding='utf-8', mode='r') as f:
				next(f)  # skip first line

				for line in f.readlines():
					####
					# SELECT *
					# FROM dataset
					# WHERE Pvalue < 0.05
					# AND Combined mutant type = ["trigenic"|"digenic"|*]

					fields = re.split(r'\t+', line)  # obtain all fields from current line separeated with tabs

					# dataset selection (now we can read from both types of dataset, s1 and s2)
					if len(fields) >= 11:
						fields.pop(5)

					# if current interaction type is not the type that we want to select, next element
					if selectedInteractionType != "*":  # will be activated if selectedInteractionType different from *
						if fields[4] != selectedInteractionType:
							continue

					if float(fields[6]) >= 0.05:  # check P-Value < 0.05, else, next
						continue

					# decide if positive or negative interaction taking in account cutoffValue
					if float(fields[5]) < cutoffValue:
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

					# Sort protoidentifier for unique key
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
				f.close()

		except ValueError as error:
			print error
		except IOError as error:
			print 'Error, file does not exist or can\'t be read'
			print error
			return 0

	# Method toString:
	#
	# Description: Returns a formatted string with all data contained in the object.
	def toString(self):

		# Returns a formatted string of a pr/npr-like 3D matrix
		def printMatrix(matrix):
			text = ''
			for i in range(self.K):
				text += str(i) + "\n\n\t\t"
				for a in range(self.K):
					text += str(a) + "\t\t\t\t\t\t\t\t"
				for j in range(self.K):
					text += str(j) + "\t\t"
					for k in range(self.K):
						for r in range(self.R):
							text += "{0:.12f}".format(matrix[i][j][k][r]) + " "
						text += "\t"
					text += "\n"
				text += "\n\n"
			return text

		# Returns a formatted string of a theta/ntheta-like two components vector
		def printVector(vector):
			text = '\n\t'
			for a in range(self.P):
				text += str(a) + "\t\t\t\t\t\t\t\t"
			for p in range(self.P):
				text += "\n"
				tab = ""
				for _ in range(2 - len(str(p)) / 4):
					tab += "\t"
				text += str(p) + tab
				for k in range(self.K):
					text += "{0:.12f}".format(vector[p][k]) + "\t"
			return text

		def printLinks(links):
			text = ''
			for r in range(self.R):
				text += "\tAparitions_R=" + str(r)
			text += '\n'
			for link in links.keys():
				tab = "\t"
				for _ in range(4 - len(link) / 4):
					tab += "\t"
				text += link + tab
				for r in range(self.R):
					text += str(links[link][r]) + "\t\t\t\t"
				text += '\n'
			return text

		text = "Likelihood: " + str(self.likelihood) + "\n"
		text += "Likelihood vector: " + str(self.vlikelihood) + "\n"
		text += "Number of genes (P): " + str(self.P) + "\n"
		text += "Number of groups of genes (K): " + str(self.K) + "\n"
		text += "Number of possible ratings (R): " + str(self.R) + "\n\n"

		# String of list of genes
		text += "LIST OF REGISTERED GENES\n"
		text += "Gene_ID\t\tGene_name\t\t\tnumAparitions\n"
		for gid in self.id_gene:
			tab, ntab = "", ""
			for _ in range(3 - (len(str(gid)) / 4)):
				tab += "\t"
			for _ in range(5 - len(self.id_gene[gid]) / 4):
				ntab += "\t"
			text += str(gid) + tab + self.id_gene[gid] + ntab + str(self.uniqueg[gid]) + '\n'

		# String of list of links by ID
		text += "\nLIST OF LINKS BETWEEN GENES ID\n"
		text += "gid1_gid2_gid3"
		text += printLinks(self.links)

		# String of list of links by gene name
		text += "\nLIST OF LINKS BETWEEN GENES ID\n"
		text += "n1_n2_n3"
		text += printLinks(self.nlinks)

		# For both pr/npr matrix
		text += "\nMATRIX OF PROBABILITIES PR\n"
		text += printMatrix(self.pr)

		text += "\nMATRIX OF PROBABILITIES NPR\n"
		text += printMatrix(self.npr)

		# Fpr both theta/ntheta vector
		text += "\n\nTHETA VECTOR\n"
		text += printVector(self.theta)

		text += "\n\nNTHETA VECTOR\n"
		text += printVector(self.ntheta)

		return text

	# Method toFile(string):
	#
	# Description: Calls the toString function and prints it in a file. Overwrites the file if the file name given is
	# identical to the file name of a file in the same folder as this script.
	#
	# Arguments:
	# 1.- Name of the output file. By-default file name will be out.txt.

	def toFile(self, filename=None):
		try:
			if filename is None:
				filename = "out.txt"
			f = codecs.open(filename, encoding='utf-8', mode="w+")
			f.write(self.toString())
			f.close()

		except IOError:
			print "I/O error"

	# Method compareLinks:
	#
	# Description: Given two object of type "Model" a and b, with a call "a.compareLinks(b)" we return a list of
	# links present in a that lack in b.
	#
	# Return parameters:
	# 1.- Returns "None" List if graph defined by "self" argument is a subgraph of "model" argument, otherwise, will
	# return the lacking links in b.

	def compareLinks(self, model):
		diff = []

		for link in self.nlinks.keys():
			if link in model.nlinks.keys():
				continue
			else:
				print link
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

	def compareGenes(self, model):
		diff = []

		for gene in self.gene_id.keys():
			if gene in model.gene_id.keys():
				continue
			else:
				diff.append(gene)
		return diff

	# Method setK(int):
	#
	# Description: Sets the K value in model (number of user groups)
	def setK(self, k):
		self.K = k

	# Method computeLikelihood:
	#
	# Description: Computes likelihood of the current data stored in the object.
	#
	# Return Parameters:
	# 1.- Returns the likelihood.
	def computeLikelihood(self):
		logL = 0.
		for triplete, ratingvector in self.links.items():  # e X ra iteration
			g1, g2, g3 = triplete.split('_')
			id1, id2, id3 = int(g1), int(g2), int(g3)  # get IDs from three genes from link
			D = [self.eps] * self.R  # Generate a vector of R position with eps value (constant defined in headers)

			for i in range(self.K):
				for j in range(self.K):
					for k in range(self.K):
						for r in range(self.R):
							D[r] += self.theta[id1][i] * self.theta[id2][j] * self.theta[id3][k] * self.pr[i][j][k][r]
			for r in range(self.R):
				logL += ratingvector[r] * math.log(D[r])
		self.likelihood = logL
		return logL

	# Method shiftValues:
	#
	# Description: Copies values from n* data structures to the current structures
	def shiftValues(self):
		self.theta = copy.copy(self.ntheta)
		for i in range(self.K):
			for j in range(self.K):
				for k in range(self.K):
					self.pr[i][j][k] = self.npr[i][j][k]

	# Method nInit:
	#
	# Description: Reinitialization of n* data structures
	def nInit(self):
		for i in range(self.P):
			self.ntheta[i] = [0.] * self.K
		for i in range(self.K):
			for j in range(self.K):
				for k in range(self.K):
					self.npr[i][j][k] = [0.] * self.R

	# Method makeIteration:
	#
	# Description: Do recursive computation to iterate over the pr, and theta vector. New values will be stored in
	# npr and ntheta instance variables. Then, new values are normalized.
	#
	# Return Parameters:
	# Updates and normalizes new values in npr and ntheta data structures.

	def makeIteration(self):

		for triplete, ratingvector in self.links.items():  # e X R iteration
			g1, g2, g3 = triplete.split('_')
			id1, id2, id3 = int(g1), int(g2), int(g3)  # get IDs from three genes from link
			D = [self.eps] * self.R  # Generate a vector of R position with eps value (constant defined in headers)

			for i in range(self.K):
				for j in range(self.K):
					for k in range(self.K):
						for r in range(self.R):
							D[r] += self.theta[id1][i] * self.theta[id2][j] * self.theta[id3][k] * self.pr[i][j][k][r]

			for i in range(self.K):
				for j in range(self.K):
					for k in range(self.K):
						for r in range(self.R):
							# auxiliary variable
							a = (self.theta[id1][i] * self.theta[id2][j] * self.theta[id3][k] * self.pr[i][j][k][r]) / D[r]
							self.ntheta[id1][i] += a * ratingvector[r]
							self.ntheta[id2][j] += a * ratingvector[r]
							self.ntheta[id3][k] += a * ratingvector[r]
							self.npr[i][j][k][r] += a * ratingvector[r]

		# Normalizations:
		# divide all possibilities of player i belonging to a group k with the number of relation of that user
		for i in range(self.P):
			for k in range(self.K):
				self.ntheta[i][k] /= float(self.uniqueg[i])

		# divide the probability of the group k giving rate to a item l with a rating r between the sum of all ratings
		for i in range(self.P):
			for j in range(self.P):
				for k in range(self.P):
					D = self.eps
					for r in range(self.R):
						D = D + self.npr[i][j][k][r]
					for r in range(self.R):
						self.npr[i][j][k][r] = self.npr[i][j][k][r] / D

	# Method compareDataset(model):
	#
	# Description: Given another object model, data of links and genes is comparated separately. This information will
	# be printed in the output.
	#
	# Return parameters:
	# 1.- True if the first dataset is subgraph for genes and links, False if not.

	def compareDataset(self, model):
		if not self.compareLinks(model):
			print "First dataset is subgraph of second dataset for links"
			node = 1
		else:
			print "First dataset is not subgraph of second dataset for links"
			node = 0
		if not self.compareGenes(model):
			print "First dataset is subgraph of second dataset for nodes"
			link = 1
		else:
			print "First dataset is not subgraph of second dataset for nodes"
			link = 0

		return link and node


# Function compareS1withS2:
#
# Description: Checks that the filtered dataset and the raw dataset are equal. The code in this function is hardcoded
# and its unique purpose is demonstrate that the selection criteria specified in article is consistent
# with obtained data from filtering dataset S1 for trigenic interaction.

def compareS1withS2():

	rawmodel = Model()
	rawmodel.getInput('Data_S1.tsv', 'trigenic', -0.08, 1)  # discards negatives for raw dataset

	treatedmodel = Model()
	treatedmodel.getInput('Data_S2.csv', 'trigenic', sys.float_info.max)  # take all data from treated dataset

	treatedmodel.toFile("treated.txt")
	rawmodel.toFile("raw.txt")
	print "\nComparing treated with raw: "
	if treatedmodel.compareDataset(rawmodel):
		sub0 = 1
	else:
		sub0 = 0

	print "\nComparing raw with treated: "
	if rawmodel.compareDataset(treatedmodel):
		sub1 = 1
	else:
		sub1 = 0

	if sub1 and sub0:
		print "\nDatasets are equal"
	else:
		print "\nDatasets are different"

# Main Function:
#
# Description: Computes the algorithm given a dataset and a number of iterations.
#
# Arguments:
# 1.- iterations: Number of iterations done by algorithm.
# 2.- samples: Number of samples done by algorithm.
# 3.- frequencyCheck: Number of iterations needed to check if likelihood has converged.
# 4.- filename: Name of the dataset filename.
# 5.- interactionType: Type of interaction selected.
# 6.- cutOffValue: Value used to determine if an interaction is positive or negative.


if __name__ == "__main__":
	random.seed(os.getpid())

	# READ ARGUMENTS
	try:
		iterations = sys.argv[1]
		samples = sys.argv[2]
		frequencyCheck = sys.argv[3]
		filename = sys.argv[4]
		interactionType = sys.argv[5]
		cutOffValue = sys.argv[6]

	# BY-DEFAULT VALUES
	except IndexError:
		iterations = 0
		samples = 10
		frequencyCheck = 25
		filename = "Data_S1.tsv"
		interactionType = "trigenic"
		cutOffValue = -0.08

	model = Model()
	model.getInput(filename, interactionType, cutOffValue, 1)

	compareS1withS2()
	print "\nFinish comparation. Beginning Iterations:"

	for sample in range(int(samples)):
		model.initializeParameters(10)
		like0 = model.computeLikelihood()
		for iteration in range(iterations):
			model.makeIteration()
			model.shiftValues()
			model.nInit()
			if iteration % frequencyCheccd Py
			k == 0:
				like = model.computeLikelihood()
				if math.fabs((like - like0) / like0) < 0.0001:
					print "Likelihood has converged"
					break
				like0 = like

		like = model.computeLikelihood()
		model.vlikelihood.append(like)  # append result into the global vector of likelihoods
		model.toFile("out"+str(sample)+".txt")

	time.sleep(5)

