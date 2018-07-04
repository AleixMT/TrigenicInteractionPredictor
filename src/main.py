#!/usr/bin/python
# -*- coding: utf-8 -*-
#################################################################################################################
# Authors: Aleix Marine i Tena (aleix.marine@estudiants.urv.cat)                                                #
# Last review: 27/6/2018                                                                                        #
# Version: 1.0                                                                                                  #
#                                                                                                               #
# Permissions: Needs permissions to read file in the same folder as the script, specified by the first argument,#
# and permissions to create a new file in the folder of the script.                                             #
#                                                                                                               #
# Parameters and descriptions: Reads the UTF-16 file in the folder where the script is located specified by the #
# first argument. Data from the first file is then written in UTF-8 format with the same name as the input file #
# plus a 'utf-8.csv' suffix. If output file exists, its content will be overwritten.                            #
#                                                                                                               #
# -Argument 1: relative route from the same folder where the script is executed.                                #
#################################################################################################################

import re
import codecs
import sys
	
from enum import Enum

class Model:

	# Constructor: Initializes data structures of algorithm

	def __init__(self):
		self.theta = []

		# BIDIRECTIONAL DICTIONARY FOR USERS
		self.id_gene = {} # dictionary that relates an id with its gene
		self.gene_id = {} # dictionary that relates a gene with its id

		# Matrix of ratings. rows: relation between gene1, gene2 and gene3 wit the format "gene1_gene2_gene3". columns: ratings (in this case 0 or 1). content: number of times seen relation between gene1, gene2 and gene3 with rating r
		self.links = {}

		# Relates the id gene with its number of interactions
		self.uniqueg = {} 
   
    # Reads data from file given by argument 1, selects type of interaction with argument 2,
    # and cutoffValue with argument 3.

	# digenic interaction cut-off in the original article (p < 0.05, |e| > 0.08)
	# trigenic interactions cut-off in the original article (p < 0.05, t < -0.08)
	# File Format:
	#
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
	# 12.- Combined mutant fitness standard deviation

	def getInput(self, filename, selectedInteractionType, cutoffValue):
		try:
			dataset = codecs.open(filename, encoding = 'utf-8', mode = 'r+')
			gid = 0
			for line in dataset:
				####
				# SELECT *
				# FROM dataset
				# WHERE Pvalue < 0.05
				# AND Combined mutant type = ["trigenic"|"digenic"]

				fields = re.split(r'\t+', line) # obtain all fields from current line
				# if current interaction type is not the type that we want to select, next element
				if fields[4] != selectedInteractionType: # first row will never be selected
					continue

				# check P-Value < 0.05, else, next element
				if float(fields[7]) >= 0.05: 
					continue
				####
				print 'khe'
				# decide if positive or negative interaction taking in account cutoffValue
				if fields[6] >= cutoffValue:
					r = 0
				else:
					r = 1

				print fields[4]
				# create a list of three elements with three genes
				print fields[1].split('+')
				print fields[3]
				gene_triplet = fields[1].split('+').extend(fields[3])
				print gene_triplet
				# REGISTER ALLELE
				for gene in gene_triplet: # for every gene in the triplet
					if gene not in self.gene_id.keys(): #  gene hasnt been seen by algorithm 
						self.gene_id[gene], n1 = gid, gid # assign a new gid to this gene 
						self.id_gene[gid] = gene # assign a new gene to this gid
						self.uniqueg[n1] = 0. # initialize this position of dictionary with "0." when user identified by n1 is the first time found
						gid += 1 # update index uid
					else: # gene is already registered
						n1 = self.gene_id[gene] # get gid from gene, already registered

					# REGISTER NUMBER OF INTERACTIONS
					self.uniqueg[n1] += 1. # indicates number of interactions for gene identified by id

				str_gene_triplet = '_'.join(gene_triplet) # joins genes with an underscore in between, indicating relation between a triplet of genes

				try:
					self.links[str_gene_triplet][r] += 1 # Indicates that link between n1 and n2 with rating r it's been seen +1 times
				except: # if link between n1 and n2 with rating r is the first time seen then
					self.links[str_gene_triplet] = [0] * 2 
					self.links[str_gene_triplet][r] += 1 #//U if tuple in links[e] has just been initialized then the operator should be assignement and not acummulation
			dataset.close()

		except:
			print 'Error, file does not exist or can\'t be read'
			return 0

	def toString(self):
		print self.id_gene
		print self.gene_id
		print self.links
		print self.uniqueg



if __name__ == "__main__":

	model = Model()
	model.getInput('aao1729_Data_S1.tsv', 'trigenic', -0.08) 
	model.toString() 
    


    

