
import sys

import _numpypy as np
import time

from math import *
#from numarray import *
#import numarray.linear_algebra as la
import copy


'''
from PyGrace.grace import Grace
from PyGrace.dataset import SYMBOLS, LINESTYLES
from PyGrace.dataset import SYMBOLS
from PyGrace.colors import ColorBrewerScheme

import sys
sys.path.append('/home/agodoy/workspace/')
sys.path.append('../../')
sys.path.append('../')

from PyGrace.grace import Grace
from PyGrace.graph import Graph
from PyGrace.dataset import DataSet
from PyGrace.Extensions.network import Network
from PyGrace.colors import ColorBrewerScheme
from PyGrace.drawing_objects import DrawText
from PyGrace.Extensions.colorbar import SolidRectangle, ColorBar
from PyGrace.axis import LINEAR_SCALE, LOGARITHMIC_SCALE
from PyGrace.drawing_objects import DrawBox
from PyGrace.drawing_objects import DrawLine
from PyGrace.Styles.metra import *
'''
import random
from math import sqrt,exp
eps= 1e-10

def maximum(vec):
	m1=0				
	for i in range(len(vec)):
		m1=max(max(vec[i]),m1)	
	return m1


def ReadData(fh,iniline=0,sep=' '):
        ##returns list of edges with number of cooperate and defect actions
        ## list of users with total number of games played
        ## list of games with number of total times the game was played
        id2user ={}
        user2id={}
        id2game={}
        game2id={}
	igot=fh.readlines()
	links={}
        uniquep = {}
        uniqueg = {}
        uid=0
        gid=0
        for line in igot[iniline:]:
             
	        about = line.strip().split(sep)
#                print about
                e1= about[0]
                e2= about[1]

                if e1 not in user2id.keys():
                        user2id[e1]=uid
                        id2user[uid]=e1
                        n1=uid
                        uniquep[n1]=0.
                        uid+=1
                else:
                        n1=user2id[e1]
                if e2 not in game2id.keys():
                        game2id[e2]=gid
                        id2game[gid]=e2
                        n2=gid
                        uniqueg[n2]=0.
                        gid+=1
                else:
                        n2=game2id[e2]
                
                uniquep[n1]+=1.
                uniqueg[n2]+=1.
                
                e= '_'.join([str(n1),str(n2)])
                ra=int(about[2])
                
                try:
                        links[e][ra]+=1
                except:
                        ratings=[0]*2
                        links[e]=ratings
                        links[e][ra]+=1
	fh.close()
        return links, uniquep,uniqueg,id2user,id2game,user2id,game2id


###########################################################3

def InitializeParameters(p,m,K,L,R):
#p - num of players; m - num of games
        
        theta=[] 
	ntheta=[] # where the new parameter values will be stored during iterations
	for i in range(p):
		a=[random.random() for _ in xrange(K)]
		theta.append(a)
		ntheta.append([0.0]*K)

	eta=[]
	neta=[]# where the new parameter values will be stored during iterations
	for i in range(m):
		a=[random.random() for _ in xrange(L)]
		eta.append(a)
		neta.append([0.0]*L)

	pr=[]
	npr=[]# where the new parameter values will be stored during iterations
	for i in range(K):
		b=[]
		c=[]
		for j in range(L):
			a=[random.random() for _ in xrange(R)]
			b.append(a)
			c.append([0.]*R)
		pr.append(b)
		npr.append(c)

	
        #Normalizations:
	for i in range(p):
		D=0.
		for k in range(K):
			D=D+theta[i][k]
		for k in range(K):
			theta[i][k]=theta[i][k]/(D+0.00000000001)

	for j in range(m):
		D=0.
		for l in range(L):
			D=D+eta[j][l]
		for l in range(L):
			eta[j][l]=eta[j][l]/(D+0.00000000001) # To�i adds this small value to avoid dividing by zero

	
	for l in range(L):
		for k in range(K):
			D=0.
			for r in range(R):
				D=D+pr[k][l][r]
			for r in range(R):
				pr[k][l][r]=pr[k][l][r]/D

        return theta,eta,pr,ntheta,neta,npr

##################################################################

def MakeIteration(theta,eta,pr,ntheta,neta,npr,K,L,R,links,users,games):

        p= len(users)
        m= len(games)
        
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

##########################################################
##########################################################
def ComputeLikelihood(links,theta,eta,pr,K,L,R):
        ##computes likelihood
        
	logL=0.
	for e,ra in links.items():
		e1,e2=e.split('_')
                n1=int(e1)
                n2=int(e2)
		D=[eps]*R
		#assuming K=L, si no hay que separar i hacer dos fors diferentes:
		for l in range(L):
				for k in range(K):
                                        for r in range(R):
                                                D[r]+=theta[n1][k]*eta[n2][l]*pr[k][l][r]

                for r in range(R):
			logL+=ra[r]*log(D[r])

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
##                   inicialització paràmetres
##                   iteraci� per actualitzar parametres 
##              

##us interessa el trainer
if __name__=="__main__":



        K=int(sys.argv[1]) #users
        L=int(sys.argv[2]) #games
        crossval=[int(a) for a in sys.argv[3:]]

        R=2	#rating 0:defect cooperate: 1  ##típicament R =5 en un recommender

        print 'K L', K,L
        random.seed(1111)
#        crossval=[1]	#5 pel 5fold cross validation
        likely=[]

        for CV in crossval:

        ## 1. llegir fitxers - no important
	        print CV
                fin = '/export/home/shared/Projects/MrBanks/gael/Data/Data MMSBM/Folded/Data_mrk+rw/DataTrain_%d_mrk+rw.csv' %(CV)
	        fh=open(fin,'r')

                links,users,games,id2u,id2g,u2id,g2id= ReadData(fh,iniline=1,sep='\t')
                p=len(users)
                m=len(games)
#                print games
	        sampling=100

	        for sam in range(sampling): ## number of times the likelihood maximizationis performed - gives better results

	                theta,eta,pr,ntheta,neta,npr= InitializeParameters(p,m,K,L,R)
			
		        Iterations=2500 ## number of iterations 
                        Like0 = ComputeLikelihood(links,theta,eta,pr,K,L,R)

		        for g in range(Iterations):
			#print g

                                ntheta,neta,npr=MakeIteration(theta,eta,pr,ntheta,neta,npr,K,L,R,links,users,games)

                                
			        theta=copy.copy(ntheta)
			        eta=copy.copy(neta)	
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
'''
#k num of parameters:
k=p*K+m*L
BIC=-2*Like+k*log(len(links))
print 'best likelihood', Like, 'num topics', -2*Like, k*log(len(links)),k, len(links)
print 'BIC?', BIC
#print 'theta', theta
'''

