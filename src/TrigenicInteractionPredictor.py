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

import re
import codecs
import os
import random
import sys
import copy
import math
#import matplotlib.pyplot as plt
import numpy as np


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

	def initialize_parameters(self, kvalue=10):
		try:
			self.K = int(kvalue)  # set K value to initialize parameters
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

			if sumtheta < self.eps: #regenerate
			        a = [random.random() for _ in range(self.K)]  # xrange to generate big lists (optimization)
			        self.theta[i] = a
                                
			sumtheta = sum(self.theta[i]) 
			for k in range(self.K):  # iterate over number of groups of genes
				try:
					self.theta[i][k] /= sumtheta
				except ZeroDivisionError:
					self.theta[i][k] /= (sumtheta + self.eps)  # adding an small value to avoid dividing by zero
#			if sum(self.theta[i]) != 1.: print ("Not 1", sum(self.theta[i]),self.theta[i])
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
	# 1.- filename: File to read will be selected with this parameters
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

	def get_input(self, argfilename, selectedinteractiontype, cutoffvalue, discard=0, numlines=0,interactions = 'ALL'):
		try:
			gid = 0
			counter = 0

			if selectedinteractiontype != 'trigenic' and selectedinteractiontype != 'digenic' and selectedinteractiontype != '*':
				raise ValueError("argument 2 selectedInteractionType must be trigenic, digenic or *")

			with codecs.open(argfilename, encoding='utf-8', mode='r') as fileref:

				line = fileref.readline()
				fields = re.split(r'\t+', line)
				if len(fields) == 12:
					reading_raw = 1
				else:
					reading_raw = 0

				for line in fileref.readlines():
					####
					# SELECT *
					# FROM dataset
					# WHERE Pvalue < 0.05
					# AND Combined mutant type = ["trigenic"|"digenic"|*]

					fields = re.split(r'\t+', line)  # obtain all fields from current line separeated with tabs

					# dataset selection (now we can read from both types of dataset, s1 and s2)
					if reading_raw:
						fields.pop(5)

					# if current interaction type is not the type that we want to select, next element
					if selectedinteractiontype != "*":  # will be activated if selectedInteractionType different from *
						if fields[4] != selectedinteractiontype:
							continue
                                                
                                                
					if interactions == 'ALL': ## decide whether we want to consider positive vs negative interactions or novel vs digenic altered interactions
                                                             
					        r = 0
					        if float(fields[6]) < 0.05 and float(fields[5]) < cutoffvalue: r = 1  # check P-Value if < 0.05, R=1 - altered interactions
						        

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

					# limit number of read lines
#					counter += 1
#					if counter == numlines:
#						break

				self.P = len(self.id_gene)  # get number of users
				fileref.close()

		except ValueError as error:
			print(error)
		except IOError as error:
			print('Error, file does not exist or can\'t be read')
			print(error)

                        
	## method that reads input links,nlinks and test_link form train and test files in the format name1_name2_name3+\t+rating  
	def get_traintest(self, trainfile,testfile): 
		try:
			gid = 0

			with codecs.open(trainfile, encoding='utf-8', mode='r') as fileref:

				for line in fileref.readlines():
					####
					# SELECT *
					# FROM dataset
					# WHERE Pvalue < 0.05
					# AND Combined mutant type = ["trigenic"|"digenic"|*]

					fields = line.strip().split('\t')  # obtain all fields from current line separeated with tabs

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

					# limit number of read lines
#					counter += 1
#					if counter == numlines:
#						break

				self.P = len(self.id_gene)  # get number of users
				fileref.close()
                                
			with codecs.open(testfile, encoding='utf-8', mode='r') as fileref:

				for line in fileref.readlines():
					####
					# SELECT *
					# FROM dataset
					# WHERE Pvalue < 0.05
					# AND Combined mutant type = ["trigenic"|"digenic"|*]

					fields = re.split(r'\t+', line)  # obtain all fields from current line separeated with tabs

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

					# Sort protoidentifier for unique key
					gene_triplet.sort()
					id_gene_triplet.sort()

					# Concatenate id and genes to create the key string
					str_name_gene_triplet = '_'.join(gene_triplet)
					str_gene_triplet = '_'.join(id_gene_triplet)  # joins genes with an underscore in between a triplet of genes

					try:
						self.test_links[str_gene_triplet][r] += 1  # link between g1, g2 and g3 with rating r it's been seen +1 times
					except KeyError:  # if link between n1 and n2 with rating r is the first time seen then
						self.test_links[str_gene_triplet] = [0] * 2
						self.test_links[str_gene_triplet][r] += 1


				self.P = len(self.id_gene)  # get number of users
				fileref.close()

		except ValueError as error:
			print(error)
		except IOError as error:
			print('Error, file does not exist or can\'t be read')
			print(error)

		print('READ DATA train',len(self.links),len(self.nlinks))
		print('READ DATA test',len(self.test_links))
					

	# Method fold:
	#
	# Description: "Folds" the data in the object model.
	# We implemented this method instead of extending the get_input method because this gives us the possibility to
	# interact with data when we know that all samples in the object are the truly selected. In the get_input method
	# we have different criteria to select samples, so probably all samples from file won't be selected.
	#
	# To do the folding we will calculate the number of samples that we have in our object model, and the we will get the
	# 20 % of the samples randomly. The samples chosen for the test set will be deleted and data structures from train_set
	# will be modified accordingly.
	#
	# Return Parameters:
	# - Fills the dictionary self.test_links with the 20 % of the links from self.links
	# - Deletes this 20 % of links from the original dictionary.
	# - Modifies the other data structures for coherence.
	def fold(self):
		test_set_size = int(len(self.links) / 5)

		# linealization of dictionary to use choice method
		arraylinks = []
		for triplet, rating in self.links.items():
			arraylinks.append(triplet)

		loop = 0 

		while loop < test_set_size:
			try:
				triplet = np.random.choice(arraylinks)  # take a triplet

				ids = triplet.split("_")  # split it

				# Obtain names for the genes in the triplet.
				names = []
				for identifier in ids:
					if self.uniqueg[int(identifier)] == 1:  # check that there's more than one aparition of that gene
						raise ValueError("Triplet "+triplet+" has at least one gene with just one aparition. Choosing randomly another")
					else:
						self.uniqueg[int(identifier)] -= 1  # substract one aparition to that gene.
						name = self.id_gene[int(identifier)]  # obtain name
						names.append(name)  # accumulate

				# obtain rating for the current triplet. (assuming just one given rating betwen r=0 or r=1)
				if self.links[triplet][0]:
					rating = 0
				else:
					rating = 1

				try:
					self.test_links[triplet][rating] += 1  # link is seen +1 time (probably the maximum will be 1)
				except KeyError:
					self.test_links[triplet] = [0] * 2  # Initialize dictionary position
					self.test_links[triplet][rating] = 1  # set 1 for the current rating

				names.sort()  # create identifier name string
				str_names = '_'.join(names)
				self.nlinks.pop(str_names)  # delete aparition from the gene name dictionary
				self.links.pop(triplet)  # delete aparition from the gene ID dictionary
				arraylinks.remove(triplet)  # delete triplet from the linearized array
				loop += 1
			except ValueError as error:
				pass

                        
	def fast_fold(self,fraction=0.2,output=None,folds=None): #if folds == 'yes' generate and print all training & test sets; fraction = relative size of test set
                
                
		test_set_size = int(len(self.links) *fraction)
		rest = len(self.links) % test_set_size
		nfolds = int(1/fraction)
                
                                

		arraylinks = []
		for triplet, rating in self.links.items():
			arraylinks.append(triplet)

                ## we shuffle the arraylinks list and then make equal splits
		np.random.shuffle(arraylinks)
                        
		if output:
			trname_file= 'train.dat'
			tname_file= 'test.dat'
			outf_train = codecs.open(trname_file, encoding='utf-8', mode="w+")
			outf_test = codecs.open(tname_file, encoding='utf-8', mode="w+")
			frange = 0 ##print only one fold
			if folds =='yes': frange = nfolds
			ta_name = []
			of_ta = []
			tra_name = []
			of_tra = []                        
			for i in range(nfolds): ##create all train test file names & open file handles
                                ta_name.append('test%d.dat' %(i))
                                of_ta.append(codecs.open(ta_name[i], encoding='utf-8', mode="w+"))
                                tra_name.append('train%d.dat' %(i))
                                of_tra.append(codecs.open(tra_name[i], encoding='utf-8', mode="w+"))                                
                        

		count = 0
		test = []
		train = []
                
		for i in range(nfolds):
			test.append(arraylinks[test_set_size*i:test_set_size*(i+1)])
			if i>0: train +=test[i]
                      
		for j in range(rest):
			print(test_set_size*nfolds+j,len(arraylinks))
			test[nfolds-1].append(arraylinks[test_set_size*nfolds+j])
                        
                         ##all remaining triplets are put inthe last test set.
                

		train += arraylinks[test_set_size*nfolds:]

		frange = 0 ##print only one fold
		if folds =='yes': frange = nfolds
                
                
                               
                ##first create & store test and train for this run
		for triplet in test[0]:
                        rating = 1
                        if self.links[triplet][0]: rating = 0
                        
                        self.test_links[triplet] = [0] * 2  # Initialize dictionary positionry:
                        self.test_links[triplet][rating] += 1  # link is seen +1 time (probably the maximum will be 1)
                        
                        ids = triplet.split("_")  # split it

			# Obtain names for the genes in the triplet.
                        names = []
                        for identifier in ids:
#				        if self.uniqueg[int(identifier)] == 1:  # check that there's more than one aparition of that gene
#					        raise ValueError("Triplet "+triplet+" has at least one gene with just one aparition. Choosing randomly another")
#				        else:
                                self.uniqueg[int(identifier)] -= 1  # substract one aparition to that gene.
                                name = self.id_gene[int(identifier)]  # obtain name
                                names.append(name)  # accumulate				

			        # obtain rating for the current triplet. (assuming just one given rating betwen r=0 or r=1)
                        names.sort()  # create identifier name string
                        str_names = '_'.join(names)
                        self.nlinks.pop(str_names)  # delete aparition from the gene name dictionary
                        self.links.pop(triplet)  # delete aparition from the gene ID dictionary


                        #print fold and train if output variable is defined
                        if output:
			
                                data = str_names + '\t'+str(rating)+'\n'
                                outf_test = of_ta[0]
                                outf_test.write(data)
                                outf_test.flush()
                                
                ##print remaining test sets
		print('remaining test sets')
                
		for f in range(1,nfolds):
		        print('fold',f,len(test[f]))
                        
		        for triplet in test[f]:
                        
                                rating = 1
                                if self.links[triplet][0]: rating = 0
                                
                                ids = triplet.split("_")  # split it


			        # Obtain names for the genes in the triplet.
                                names = []
                                for identifier in ids:

                                        name = self.id_gene[int(identifier)]  # obtain name
                                        names.append(name)  # accumulate				

			        # obtain rating for the current triplet. (assuming just one given rating betwen r=0 or r=1)
                                names.sort()  # create identifier name string
                                str_names = '_'.join(names)

                                #print fold and train if output variable is defined
                                if output:
			
                                        data = str_names + '\t'+str(rating)+'\n'
                                        outf_test = of_ta[f]
                                        outf_test.write(data)
                                        outf_test.flush()


		for triplet in train:
			count +=1 
			if output:
                                ids = triplet.split("_")  # split it
				        # Obtain names for the genes in the triplet.
                                names = []
                                for identifier in ids:
                                        name = self.id_gene[int(identifier)]
                                        names.append(name)
					
                                                
				        # obtain rating for the current triplet. (assuming just one given rating betwen r=0 or r=1)
                                rating = 1
                                if self.links[triplet][0]:
                                        rating = 0
						
                                names.sort()
                                str_names = '_'.join(names)
                                data = str_names + '\t'+str(rating)+'\n'
                                outf_train = of_tra[0]
                                outf_train.write(data)
                                outf_train.flush()
					
		outf_train.close()
		outf_test.close()
		

                                        
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
#			print(triplet,rating,self.test_links[triplet])
                        
		# Sort results from less to higher predicted probability
		self.results.sort()
		self.results.reverse()

	# Method calculate_metrix:
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

		# Prints the likelihood vector in csv format splitted by tabs (\t)
		def tostring_likelihood(vector):
			txt = "Sample\titeration\tlikelihood\n"
			for num_sample, num_iteration, num_likelihood in vector:
				txt += str(num_sample) + "\t" + str(num_iteration) + "\t" + str(num_likelihood) + "\n"
			return txt

		# Returns a csv string of a pr/npr-like 3D matrix
		def print_matrix(matrix):
			txt = ''
			for i in range(self.K):
				txt += str(i) + "\n\t"
				for a in range(self.K):
					txt += str(a) + "\t\t"
				txt += "\n\t"
				for b in range(self.K):
					txt += "0\t1\t"
				txt += "\n"
				for j in range(self.K):
					txt += str(j) + "\t"
					for k in range(self.K):
						for r in range(self.R):
							txt += "{0:.6f}".format(matrix[i][j][k][r]) + "\t"
					txt += "\n"
				txt += "\n\n"
			return txt

		# Returns a formatted string of a theta/ntheta-like two components vector
		def print_vector(vector):
			txt = '\n\t'
			for a in range(self.K):
				txt += str(a) + "\t"
			for p in range(self.P):
				txt += "\n"
				txt += str(p) + "\t"
				for k in range(self.K):
					txt += "{0:.12f}".format(vector[p][k]) + "\t"
			txt += '\n'
			return txt

		def print_links(links):
			txt = ''
			for r in range(self.R):
				txt += "\tAparitions_R=" + str(r)
			txt += '\n'
			for link in links.keys():
				txt += link + "\t"
				for r in range(self.R):
					txt += str(links[link][r]) + "\t"
				txt += '\n'
			return txt

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
		text += "Likelihood vector: \n" + str(tostring_likelihood(self.vlikelihood))
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

		# String of list of links by ID
		text += "\nLIST OF LINKS BETWEEN GENE IDS\n"
		text += "gid1_gid2_gid3"
		text += print_links(self.links)

		# String of list of links by gene name
		text += "\nLIST OF LINKS BETWEEN GENE NAMES\n"
		text += "n1_n2_n3"
		text += print_links(self.nlinks)

		# For both pr/npr matrix
		text += "\nMATRIX OF PROBABILITIES PR\n"
		text += print_matrix(self.pr)

		text += "\nMATRIX OF PROBABILITIES NPR\n"
		text += print_matrix(self.npr)

		# Fpr both theta/ntheta vector
		text += "\nTHETA VECTOR\n"
		text += print_vector(self.theta)

		text += "\nNTHETA VECTOR\n"
		text += print_vector(self.ntheta)

		return text
	# Method tostring:
	#
	# Description: Returns a CSV-like (using one tab as separator between fields) format string with data from the model
	# object.
	def to_string_short(self):

		# Prints the likelihood vector in csv format splitted by tabs (\t)
		def tostring_likelihood(vector):
			txt = "Sample\titeration\tlikelihood\n"
			for num_sample, num_iteration, num_likelihood in vector:
				txt += str(num_sample) + "\t" + str(num_iteration) + "\t" + str(num_likelihood) + "\n"
			return txt

		# Returns a csv string of a pr/npr-like 3D matrix
		def print_matrix(matrix):
			txt = ''
			for i in range(self.K):
				txt += str(i) + "\n\t"
				for a in range(self.K):
					txt += str(a) + "\t\t"
				txt += "\n\t"
				for b in range(self.K):
					txt += "0\t1\t"
				txt += "\n"
				for j in range(self.K):
					txt += str(j) + "\t"
					for k in range(self.K):
						for r in range(self.R):
							txt += "{0:.6f}".format(matrix[i][j][k][r]) + "\t"
					txt += "\n"
				txt += "\n\n"
			return txt

		# Returns a formatted string of a theta/ntheta-like two components vector
		def print_vector(vector):
			txt = '\n\t'
			for a in range(self.K):
				txt += str(a) + "\t"
			for p in range(self.P):
				txt += "\n"
				txt += str(p) + "\t"
				for k in range(self.K):
					txt += "{0:.12f}".format(vector[p][k]) + "\t"
			txt += '\n'
			return txt

		def print_links(links):
			txt = ''
			for r in range(self.R):
				txt += "\tAparitions_R=" + str(r)
			txt += '\n'
			for link in links.keys():
				txt += link + "\t"
				for r in range(self.R):
					txt += str(links[link][r]) + "\t"
				txt += '\n'
			return txt

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
#		text += "Likelihood vector: \n" + str(tostring_likelihood(self.vlikelihood))
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

	def to_file_short(self, name_file=None):
		try:
			if name_file is None:
				name_file = "out.txt"
			fileref = codecs.open(name_file, encoding='utf-8', mode="w+")
			data = self.to_string_short()
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
	# Description: Computes likelihood of the current data stored in the object.
	#
	# Return Parameters:
	# 1.- Returns the likelihood.
	def compute_likelihood(self):
		log_l = 0.
		for triplet, rating_vector in self.links.items():  # e X ra iteration
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
                        
#			s1 = sum(self.theta[id1])
#			s2 = sum(self.theta[id2])
#			s3 = sum(self.theta[id3])
			counter[id1]+=1
			counter[id2]+=1
			counter[id3]+=1

#			txt = "Sums theta %f %f %f" % (s1,s2,s3)

                        
#			print (txt)

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
##				self.ntheta[i][k] /= float(self.uniqueg[i]) ## wrong denominator, uniqueg was computed for the whole set!! Needs to be computed only for the train set
				self.ntheta[i][k] /= float(counter[i])
##			s1 = sum(self.ntheta[i])
##			s2 = sum(self.theta[id2])a
##			s3 = sum(self.theta[id3])
##			txt = "Sums theta %f %f %f" % (s1,s2,s3)
##			txt = "Sum ntheta %f " % (s1)
##			if self.uniqueg[i]!= counter[i]: print ("Counter not equal",self.uniqueg[i],counter[i],i)
##			print (txt)
                         
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

	# Method plot_likelihood:
	#
	# Description: Plots the likelihood.
	def plot_likelihood(self):
		i = 0
		print(self.vlikelihood[0][0])
		print(self.vlikelihood)
		data_x = [self.vlikelihood[0][1]]
		data_y = [self.vlikelihood[0][2]]
		print(self.vlikelihood[0][2])

		itera = self.vlikelihood[0][1]
		while itera < self.vlikelihood[1+i][0]:
			data_x.append(self.vlikelihood[1+i][1])
			print("bucle x")
			print(data_x)
			data_y.append(self.vlikelihood[1+i][2])
			print("bucle y")
			print(data_y)

		plt.plot(data_x, data_y)
		plt.title('likelihood over iterations')
		plt.show()

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
		iterations = int(sys.argv[1])
		samples = int(sys.argv[2])
		frequencyCheck = int(sys.argv[3])
		filename = sys.argv[4]
		interactionType = sys.argv[5]
		cutOffValue = float(sys.argv[6])
		argk = int(sys.argv[7])

	# BY-DEFAULT VALUES
	except IndexError:
		iterations = 10000
		samples = 100
		frequencyCheck = 10
		filename = "/export/home/shared/Projects/Trigenic/TrigenicInteractionPredictor/src/Data_S1.csv"
		interactionType = "trigenic"
		cutOffValue = -0.08
		argk = 2
	
	msg = "\n****************************************\n* Trigenic Interaction Predictor v 1.0 *\n**************"
	msg += "**************************\n\nDoing "+str(samples)+" samples of "+str(iterations)+" iterations."
	msg += "\nData is read from file "+filename+"."+"\n"+interactionType+" interactions are currently selected. "
	msg += "\nTau/epsilon cutOffvalue is "+str(cutOffValue)+"K value (number of groups) is "+str(argk)+"."
	msg += "\nLikelihood will be calculated every "+str(frequencyCheck)+" iterations."
	print(msg)

	model = Model()
	model.get_input(filename, interactionType, cutOffValue, 0, 10000)
	model.fast_fold(output=1,folds='yes')

	print("\nStarting algorithm...")

	for sample in range(int(samples)):
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
