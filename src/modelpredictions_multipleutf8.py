import sys, os
import time
from math import *
import copy
import random
from math import sqrt,exp

###
# DICTIONARY MARKS:
# /// --> TO-DO, get to it later, needs revision...
# //RF --> Re-factored
# //rf --> may need refactoring
# //R --> reference
# //M --> Modified from original code
# //? --> unknown behavior
# //U --> unexpected behaviour or statement
# //RTFM --> Read The Fuckin Manual
# //* --> Start dictionary mark
# //** --> End dictionary mark
###

eps= 1e-10 #//?

# Returns maximum of a given vector
def maximum(vec):
	m1=0				
	for i in range(len(vec)):
		m1=max(max(vec[i]),m1)	
	return m1

## returns list of edges with number of cooperate and defect actions
## list of users with total number of games played
## list of games with number of total times the game was played

# sep is argument for separator in file
# iniline is the first line readed
# fh is reference to file

# //rf iniline is an argument not a parameter
# //rf sep is an arguemnt not a parameter. sep = ' ' is the default behaviour. Taking in account low reusability is not necessary //rf

def ReadData(fh,iniline=0,sep=' '):
        # BIDIRECTIONAL DICTIONARY FOR USERS
        id2user = {} # dictionary that relates user with id
        user2id = {} # dictionary that relates id with user

        # BIDIRECTIONAL DICTIONARY FOR GAMES
        id2game = {} # dictionary that relates game with id
        game2id = {} # dictionary that relates id with games
        
        igot = fh.readlines() #//rf content from dataset
        
        links = {} # Matrix of ratings. rows: relation between u and g with the format "uid_gid". columns: ratings (in this case 0 or 1). content: number of times seen relation between u and g with rating r #//U//R on line 99
        
        uniquep = {} # Relates the id user with its number of interactions
        uniqueg = {} # Relates the id user with its number of interactions
        
        uid = 0 # autocorrelative number given to a user
        gid = 0 # autocorrelative number given to a game
        
        for line in igot[iniline:]: # read file line by line
             
                u, g, r = line.strip().split(sep) #//M//RF Variable name is not adequate # Obtain user u, game g and rating r. 

                #//rf* inheritance with dynamic binding method or wrapper class
                # REGISTER USER
                if u not in user2id.keys(): # user u hasnt been seen by algorithm #//rf "not in" preceded by "else" statements should be written as "in" "else" for optimization
                        user2id[u] = uid #//M assign a new uid to this user 
                        id2user[uid] = u #//M assign a new user to this uid
                        n1 = uid #//rf temporary variable
                        uniquep[n1] = 0. # initialize this position of dictionary with "0." when user identified by n1 is the first time found
                        uid += 1 # update index uid
                else: # user u is already registered
                        n1 = user2id[u] #//M get uid for user u, already registered

                # REGISTER GAME
                if g not in game2id.keys():
                        game2id[g] = gid #//M
                        id2game[gid] = g #//M
                        n2 = gid
                        uniqueg[n2] = 0.
                        gid += 1
                else:
                        n2 = game2id[g] #//M

                #//rf**
                
                #// REGISTER NUMBER OF INTERACTIONS
                uniquep[n1] += 1. # indicates number of interactions for user identified by n1 id
                uniqueg[n2] += 1.
                
                e = '_'.join([str(n1), str(n2)]) # joins n1 and n2 with an underscore in between, indicating relation between user identified by n1 and game identified by n2
                
               # r = int(r) #//M//RF 
                
                try:
                        links[e][int(r)] += 1 #//M Indicates that link between n1 and n2 with rating r it's been seen +1 times
                except: # if link between n1 and n2 with rating r is the first time seen then
                        links[e] = [0] * 2 # //RF ratings is useless variable since statement has pipe behaviour # //? if r is a rating, then link[e] should be a tuple of 5 elements #//R line 334. Model behavies as a like/doesn't like ratings
                        links[e][int(r)] += 1 #//U if tuple in links[e] has just been initialized then the operator should be assignement and not acummulation
        fh.close()
        return links, uniquep, uniqueg, id2user, id2game, user2id, game2id


###########################################################3
# Initializes theta, eta, pr, ntheta, neta, npr
# theta --> vector of possibilities of a user belonging to a determinate group of users
# eta --> vector of possibilities of a user belonging to a determinate group of users
# Given p, m, K, L, R values
# p --> num of players
# m --> num of games
# K --> number of groups of players
# L --> number of groups of items
# R --> Number of possible ratings
def InitializeParameters(p, m, K, L, R):
        
    theta = [] 
    ntheta = [] # where the new parameter values will be stored during iterations
    for i in range(p): # Iterate over the number of players
        a = [random.random() for _ in xrange(K)] # xrange to generate big lists
        theta.append(a) # appends to theta a vector of random values
        ntheta.append([0.0] * K) # generate a vector of reals with number of players size and append it to ntheta 

    eta = []
    neta = [] # where the new parameter values will be stored during iterations
    for i in range(m): # Iterate over the number of games
        a = [random.random() for _ in xrange(L)]
        eta.append(a) # appends to eta a vector of random values
        neta.append([0.0] * L) # generate a vector of reals with number of games size and append it to ntheta 

    pr = []
    npr = [] # where the new parameter values will be stored during iterations
    for i in range(K): # Iterate over the number of groups of players
        b = []
        c = []
        for j in range(L): # iterate over the number of groups of items
            a = [random.random() for _ in xrange(R)] # Generate a vector of random values with the size of the number of possible ratings
            pr.append(b.append(a)) # //M*//RF* pipe-like structure //U pr is a repetitive diagonal matrix [[a[1], [a[2]], [a[3]], ... , [a[L]]] [K] Is this probably a mistake?
            npr.append(c.append([0.] * R)) # //** //U npr is a repetitive diagonal matrix [[0.[R][1]], [0.[R][2]], [0.[R][3]], ... , [0.[R][L]]] [K]
    # //rf* duplicated code
    # Normalization for theta vector of players:
    for i in range(p): # Iterate over number of players
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
                D = D + pr[k][l][r] # //R line 138. We access pr matrix as if it was squared (I mean, complete. Probably is a mistake 
            for r in range(R):
                pr[k][l][r] = pr[k][l][r] / D # sum of probabilities of a user from group K giving rate R to a item from group L is 1

        return theta, eta, pr, ntheta, neta, npr

##################################################################

def MakeIteration(theta, eta, pr, ntheta, neta, npr, K, L, R, links, users, games):

    p = len(users)
    m = len(games)
        
    for e, ra in links.items(): # e X ra iteration
        e1, e2 = e.split('_')
        n1 = int(e1)
        n2 = int(e2) # get identifiers of item and user related
        D = [eps] * R
        for l in range(L):
            for k in range(K):
                   for r in range(R):
                        D[r] += theta[n1][k] * eta[n2][l] * pr[k][l][r] ## auxiliry variable normalization denominator of omega in paper
                
        for l in range(L):
            for k in range(K):
                   for r in range(R):
                        a=(theta[n1][k]*eta[n2][l]*pr[k][l][r])/D[r] ## auxiliary variable
                        ntheta[n1][k]+=a*ra[r]
                        neta[n2][l]+=a*ra[r]
                        npr[k][l][r]+=a*ra[r]
    #Normalizations:
    for i in range(p):
        for k in range(K):
            ntheta[i][k]/=float(users[i])

    for j in range(m):
        for l in range(L):
            neta[j][l]/=float(games[j])

    for l in range(L):
        for k in range(K):
            D=eps
            for r in range(R):
                D=D+npr[k][l][r]
            for r in range(R):
                npr[k][l][r]=npr[k][l][r]/D
                
        return ntheta,neta,npr


##########################################################
def MakeIterationPrior(theta,eta,pr,ntheta,neta,npr,K,L,R,links,alpha):

        corr={}

	for j in range(m):
		corr[j]=[0.]*L
		for nei in games[j]:
			F=0.
			for l in range(L):
				F+=eta[j][l]*eta[nei][l]
				Ftot+=(1.-F)
				#print (1.-F)
				for l in range(L):
					corr[j][l]+=F-eta[nei][l]
        
        for e,ra in links.items():
		e1,e2=e.split('_')
                n1=int(e1)
                n2=int(e2)
		D=[eps]*R
#                print theta[n1]
		for l in range(L):
			for k in range(K):
                                for r in range(R):
					D[r]+=theta[n1][k]*eta[n2][l]*pr[k][l][r] ## auxiliry variablenormalization
                
		for l in range(L):
			for k in range(K):
                                for r in range(R):
					a=(theta[n1][k]*eta[n2][l]*pr[k][l][r])/D[r] ## auxiliary variable

					ntheta[n1][k]+=a*ra[r]
					neta[n2][l]+=a*ra[r]
					npr[k][l][r]+=a*ra[r]
			
		#Normalizations:

	for i in range(p):
                for k in range(K):
			ntheta[i][k]/=float(users[i])

        for j in range(m):

		for l in range(L):
			neta[j][l]/=float(games[j])
			
	for l in range(L):
		for k in range(K):
			D=eps
			for r in range(R):
				D=D+npr[k][l][r]
			for r in range(R):
				npr[k][l][r]=npr[k][l][r]/D
#	print n1, ntheta[n1]	

        return ntheta,neta,npr

##################################################################################
# Likelihood defines how good our model and parameters describe our data.
# is defined as 
# for every u rating and item i from the observed data
#     sum of logs of
#          sum over k groups of users and l groups of items and r possible ratings
#                theta[user][k] * eta[item][l] * pr[k][l][r]
##################################################################################
def ComputeLikelihood(links, theta, eta, pr, K, L, R):
        ##computes likelihood
        
    logL = 0.
    for e, ra in links.items(): # e X ra iteration
        e1, e2 = e.split('_')
        n1, n2 = int(e1), int(e2) # //RF//M get IDs (player and item) from link
        D = [eps] * R # Generate a vector of R position with eps value (constant defined in headers)
        
        #assuming K=L, si no hay que separar i hacer dos fors diferentes:
        for l in range(L):
            for k in range(K):
                for r in range(R):
                    D[r] += theta[n1][k] * eta[n2][l] * pr[k][l][r]
        for r in range(R):
            logL += ra[r] * log(D[r])

        return logL

###############################################################

def PrintResults(fname,logL,theta,eta,pr,nusers,ngames,K,L,R):
        
	fout=open(fname,'a')
	fout.write('%s\n' % logL)

        for i in range(nusers):
		for k in range(K):
			fout.write('%s ' % theta[i][k])
		fout.write('\n')
	for j in range(ngames):
		for l in range(L):
			fout.write('%s ' % eta[j][l])
		fout.write('\n')
	for k in range(K):
		for l in range(L):
                        for r in range(R):
			        fout.write('%s ' % (pr[k][l][r]))
		        fout.write('\n')
	fout.close()

        return
##################################################################


## pseudo code del qe ve ara:
## bucle per tots els crossvalidations
##     bucle per mostres
##              Trainer
##                   inicialitzaciÃ³ parÃ metres
##                   iteració per actualitzar parametres 
##              

##us interessa el trainer
if __name__=="__main__":

        random.seed(os.getpid()) #//M//RF Fixed-value seed initialization does not make sense
        
        K = int(sys.argv[1]) # Number of groups of users
        L = int(sys.argv[2]) # Number of groups of games
        R = 2 #rating: 0 --> defect, 1 --> cooperate. Typically R = 5 in a recommender 
        
        likely = []
        crossval = [int(a) for a in sys.argv[3:]] # argument 3 or more are converted into int, that are crossvalidation values

        for CV in crossval:
                # links[p X m][R]:  Matrix of ratings. rows: relation between u and g with the format "uid_gid". columns: ratings (in this case 0 or 1). content: number of times seen relation between u and g with rating r #//U//R on line 99
                # users: Relates the id user with its number of interactions
                # games: Relates the id game with its number of interactions
                # id2u: Dictionary that relates user with id
                # id2g: Dictionary that relates game with id
                # u2id: Dictionary that relates id with user
                # g2id: Dictionary that relates id with games
                links, users, games, id2u, id2g, u2id, g2id = ReadData(open('/export/home/shared/Projects/MrBanks/gael/Data/Data MMSBM/Folded/Data_mrk+rw/DataTrain_%d_mrk+rw.csv' % (CV), 'r'), iniline = 1, sep = '\t') # //RF generate filename of the datatrain used, create proxy to desired file and get data from it
                
                p = len(users) # p: number of users
                m = len(games) # m: number of games
                sampling = 100 # number of times the likelihood maximizationis performed - gives better results

                for sam in range(sampling):
                    # theta: theta[p][K]. matrix that relates a player p belonging to group K with its probability of belonging to it
                    # eta: eta[m][L] matrix that relates an item m belonging to group L with its probability of belonging to it 
                    # pr: pr[k][l][r] 3D matrix that relates a group of users k, group of items l with the probability with the probability of giving a rate r
                    # ntheta, neta and npr have the same dimensions as its counterparts, but all of his values are initializated to zero
                    theta, eta, pr, ntheta, neta, npr = InitializeParameters(p, m, K, L, R)
                    Iterations = 2500 # number of iterations 
                    Like0 = ComputeLikelihood(links, theta, eta, pr, K, L, R) # Obtain Likelihood of our random initializated values

                    for g in range(Iterations): # make Iterations iterations
                        ntheta, neta, npr = MakeIteration(theta, eta, pr, ntheta, neta, npr, K, L, R, links, users, games)
                        theta = copy.copy(ntheta)
                        eta = copy.copy(neta)
                        for k in range(K):
                            for l in range(L):
                                pr[k][l]=npr[k][l]
                        for i in range(p):
                            ntheta[i]=[0.]*K
                        for j in range(m):
                            neta[j]=[0.]*L
                        for k in range(K):
                            for l in range(L):
                                npr[k][l]=[0.]*R
                        if g % 25 ==0:
                            Like = ComputeLikelihood(links,theta,\
                                                     eta,pr,K,L,R)
                            if fabs((Like-Like0)/Like0)<0.0001: break
                            Like0=Like

                                        
                    Like = ComputeLikelihood(links,theta,eta,pr,K,L,R)
                    likely.append(Like)
                    fname='testCV%dK%dL%d.dat' %(CV,K,L)
                    PrintResults(fname,Like,theta,eta,pr,\
                                     len(users),len(games),K,L,R)

                    print sam, Like,g

                print 'best likelyhood', likely.index(max(likely)),max(likely),len(users),len(games),K,L,R

                fout= 'userid%d.dat' %(CV)
                fo=open(fout,'w')
                for id,u in id2u.items():
                        print>> fo,id,u
                
                fo.close()
                fout= 'gameid%d.dat' %(CV)
                fo=open(fout,'w')
                for id,g in id2g.items():
                        print>> fo,id,g

                fo.close()
