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
import io
import os
import random	

class Model:

	# Constructor: Initializes data structures of algorithm

	def __init__(self):
                # Vectors of probability, that relates the probability of user u belonging to a determinate group
                self.ntheta = []
		self.theta = []

		# BIDIRECTIONAL DICTIONARY FOR USERS
		self.id_gene = {} # dictionary that relates an id with its gene
		self.gene_id = {} # dictionary that relates a gene with its id

                # Probability matrix
                self.pr = []
                self.npr = []
                
		# Matrix of ratings. rows: relation between gene1, gene2 and gene3 wit the format "gene1_gene2_gene3". columns: ratings (in this case 0 or 1). content: number of times seen relation between gene1, gene2 and gene3 with rating r
		self.links = {}

		# Relates the id gene with its number of interactions
		self.uniqueg = {} 

                # R: number of possible ratings (tags) in a link. Typically 5 in a recommender. In our case, 1 for interaction 0 for no interaction
                self.R = 2

                # K: Number of groups
                self.K = 0

                # P: Number of genes
                self.P = 0
   
        ###########################################################3
        # Initializes theta, pr, ntheta, npr
        # theta --> vector of possibilities of a gene belonging to a determinate group of genes
        # Given P, K, R values
        # P --> num of genes
        # K --> number of groups of genes
        # R --> Number of possible ratings
        def InitializeParameters(self):
                for i in range(self.P): # Iterate over the number of players
                        a = [random.random() for _ in xrange(K)] # xrange to generate big lists (optimization)
                        self.theta.append(a) # appends to theta a vector of random values
                        self.ntheta.append([0.0] * self.K) # generate a vector of reals with number of players size and append it to ntheta 

                for i in range(self.K):
                        b = []
                        c = []
                        for j in range(self.K):
                                d = []
                                e = []
                                for k in range(self.K):
                                        d.append(random.random() for _ in xrange(R))
                                        e.append([0.] * R)
                                b.append(d)
                                c.append(e)
                        self.pr.append(b)
                        self.npr.append(c)
                                        
                # Normalization for theta vector of players:
                for i in range(self.P): # Iterate over number of players
        D = 0. 
        for k in range(K): # iterate over number of groups of players
            D = D + theta[i][k] 
        # D = sum of possibilities of theta vector for user i
        for k in range(K): 
            theta[i][k] = theta[i][k] / (D + 0.00000000001)

    # Normalization for eta vector of items:
    for j in range(m):
        D = 0.
        for l in range(L):
            D = D + eta[j][l]
        for l in range(L):
            eta[j][l] = eta[j][l] / (D + 0.00000000001) # Toñi adds this small value to avoid dividing by zero
    # //**
    for l in range(L):
        for k in range(K):
        D = 0.
            for r in range(R):
                D = D + pr[k][l][r] 
            for r in range(R):
                pr[k][l][r] = pr[k][l][r] / D # sum of probabilities of a user from group K giving rate R to a item from group L is 1

        return theta, eta, pr, ntheta, neta, npr


        # File Format:
	#
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
    # digenic interaction cut-off in the original article (p < 0.05, |e| > 0.08)
    # trigenic interactions cut-off in the original article (p < 0.05, t < -0.08)
	# 1.- Query strain ID
	# 2.- Query allele name
	# 3.- Array strain ID
	# 4.- Array allele name
	# 5.- Combined mutant type
        # 6.- Adjusted genetic interaction
        # 7.- P-value
        # 8.- Interaction type
    #
    # Method getInput:
    #
    # Reads data from file in the same folder as this script with name given by argument 1 "filename"
    # selects type of interaction with argument 2 "selectedInteractionType" [trigenic|digenic|*]
    # and selects interaction as positive [1] with an adjusted genetic interaction score under "cutoffValue" given by argument 3.
    # selected as negative [0] if not. Use sys.float_info.max or sys.float_info.min for assuming all positive or negative interaction
    # All genes will be selected with P-value < 0.05
    # Specify which dataset will be read with argument 4 "dataset" [1|2]


    
        def getInput(self, filename, selectedInteractionType, cutoffValue, dataset):
	        try:
			gid = 0
            
                        if int(dataset) != 1 and int(dataset) != 2:
                                raise ValueError("argument 4 dataset must be 1 or 2")

                        if selectedInteractionType != 'trigenic' and selectedInteractionType != 'digenic' and selectedInteractionType != '*':
                                raise ValueError("argument 2 selectedInteractionType must be trigenic, digenic or *")

			with codecs.open(filename, encoding = 'utf-8', mode = 'r') as file:
                                next(file) # skip first line
                
                                for line in file.readlines():
                                        ####
			                # SELECT *
			                # FROM dataset
			                # WHERE Pvalue < 0.05
			                # AND Combined mutant type = ["trigenic"|"digenic"|*]

                                        fields = re.split(r'\t+', line) # obtain all fields from current line
                                
                                        # dataset selection (now we can read from both types of dataset, s1 and s2)
                                        if dataset == 1:
                                                fields.pop(5)
                                        
                                        # if current interaction type is not the type that we want to select, next element
                                        if selectedInteractionType != "*": # selection condition will be activated if selectedInteractionType different from *
                                                if fields[4] != selectedInteractionType: 
                                                        continue
                                
                                        if float(fields[6]) >= 0.05: # check P-Value < 0.05, else, next else
                                                continue

                
                                        # decide if positive or negative interaction taking in account cutoffValue
                                        if float(fields[5]) < cutoffValue:
                                                r = 1
                                        else:
                                                r = 0

                                        # create list with the three implicated alleles
                                        gene_triplet = fields[1].split('+') 
                                        gene_triplet.append(fields[3])
                                        gene_triplet.sort() # sort genes in order to have an unique criteria and recognize repetitions

                                        # REGISTER ALLELE
                                        for gene in gene_triplet: # for every gene in the triplet
                                                if gene not in self.gene_id.keys(): #  gene hasnt been seen by algorithm
                                                        self.gene_id[gene], n1 = gid, gid # assign a new gid to this gene 
                                                        self.id_gene[gid] = gene # assign a new gene to this gid
                                                        self.uniqueg[n1] = 0 # initialize this position of dictionary with "0." when user identified by n1 is the first time found //U should be an array of integers
                                                        gid += 1 # update index uid
                                                else: # gene is already registered
                                                        n1 = self.gene_id[gene] # get gid from gene, already registered

                                                # REGISTER NUMBER OF INTERACTIONS
                                                self.uniqueg[n1] += 1 # indicates number of interactions for gene identified by id //U should be integers

                                        str_gene_triplet = '_'.join(gene_triplet) # joins genes with an underscore in between, indicating relation between a triplet of genes

			                try:
				                self.links[str_gene_triplet][r] += 1 # Indicates that link between n1 and n2 with rating r it's been seen +1 times
		                        except: # if link between n1 and n2 with rating r is the first time seen then
                                                self.links[str_gene_triplet] = [0] * 2 
			                        self.links[str_gene_triplet][r] += 1 #//U if tuple in links[e] has just been initialized then the operator should be assignement and not acummulation
                        self.P = len(self.id_gene)
                        file.close()
                except ValueError as e:
                        print e
		except IOError as e:
			print 'Error, file does not exist or can\'t be read'
			return 0

	def toString(self):
                print "LIST OF REGISTERED GENES"
                for id in self.id_gene:
                        print "Gene with id "+str(id)+' is '+self.id_gene[id]
                        print "And has "+str(self.uniqueg[id])+" aparitions in dataset."
                        print

                print "LIST OF LINKS BETWEEN GENES"
                for link in self.links:
                        print
                        print "Triplete of genes "+link
                        for r in range(0, self.R):
                                print "Has "+str(self.links[link][r])+" aparitions with rate "+str(r)+"."

        def toFile(self):
                try:
                        f = codecs.open("genes.txt", encoding = 'utf-8', mode = "w+")
                        f.write("Gene_ID Gene_name numAparitions\n")
                        for id in self.id_gene:
                                f.write (str(id)+" "+self.id_gene[id]+" "+str(self.uniqueg[id])+"\n")

                        f.close()
                        
                        f = codecs.open("links.txt", encoding = 'utf-8', mode = "w+")
                        f.write("Triplete_Name ")
                        for r in range (0, self.R):
                                f.write("rating="+str(r)+" ")
                        f.write("\n")
                        
                        for link in self.links:
                                f.write(link)
                                for r in range (0, self.R):
                                        f.write(" "+str(self.links[link][r]))
                                f.write("\n")
                        f.close()
                        
                except IOError as e:
                        print "I/O error"

        # Method compareLinks
        # Given two object of type "Model" a and b, with a call "a.compareLinks(b)" we return a list of links present in a that lack in b.
        # --> Returns "None" List if graph defined by "self" argument is a subgraph of "model" argument
        def compareLinks(self, model):
                diff = []
                
                for link in self.links:
                        if link in model.links:
                                pass
                        else:
                                diff.append(link)
                return diff

        # Method compareGenes
        # Given two object of type "Model" a and b, with a call "a.compareLinksGenes(b)" we return a list of genes present in a that lack in b.
        # --> Returns "None" List if graph defined by "self" argument is a subgraph of "model" argument
        def compareGenes(self, model):
                diff = []

                for gene in self.gene_id.keys():
                        if gene in model.gene_id.keys():
                                pass
                        else:
                                diff.append(gene)
                return diff                

        # Method setK
        # Sets the K value in model (number of user groups)
        def setK(self, K):
                self.K = K
        
                
                
def compareS1withS2():
        rawmodel = Model()
	rawmodel.getInput('Data_S1.tsv', 'trigenic', -0.08, 1)
        treatedmodel = Model()
        treatedmodel.getInput('Data_S2.csv', 'trigenic', sys.float_info.max, 2)
        print treatedmodel.compareLinks(rawmodel)
        print treatedmodel.compareGenes(rawmodel)
        
        print "S1 trigenic interactions with t < -0.08 and P-value < 0.05 "
        if not treatedmodel.compareLinks(rawmodel):
                print "have equal links "
        else:
                print "have different links "
        if not treatedmodel.compareGenes(rawmodel):
                print "have equal genes"
        else:
                print "have different genes"
        print "as trigenic interaction in S2 dataset"
        
        
if __name__ == "__main__":
        random.seed(os.getpid())
        
        model = Model()
        model.getInput("Data_S1.tsv", "trigenic", -0.08, 1)

        model.setK(int(sys.argv[1]))
    


    

