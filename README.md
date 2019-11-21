# Trigenic Interaction Predictor
Biotechnology Final Degree Project which consists of an algorithm that predicts the probability of interaction between triplets of genes from model organism *[Pichia pastoris](https://en.wikipedia.org/wiki/Pichia_pastoris)* (*P. pastoris*) using Mixed-Membership [Stochastic Block Model](https://en.wikipedia.org/wiki/Stochastic_block_model) (MMSBM).

## Features
* We use supplementary materials from the article [Systematic analysis of complex genetic interactions"](http://science.sciencemag.org/content/360/6386/eaao1729) DOI: 10.1126/science.aao1729 to get our datasets. These datasets can also be found in the `doc/` folder in this repository.
* Due to the use of big datasets in our project we included [git-lfs](https://git-lfs.github.com/) in order to work with big files comfortably.
* The main code is written is Python3 and can be run using the standard interpreter, but for optimization purposes we recommend to use [pypy3](https://pypy.org/) for [Ubuntu Linux x64](https://bitbucket.org/pypy/pypy/downloads/pypy3.6-v7.1.1-linux64.tar.bz2).
* The training algorithm can be executed in parallel using [GNU-parallel](https://www.gnu.org/software/parallel/).
* To write the corresponding article we used [LaTeX](https://www.latex-project.org/).

## Introduction
Genetic interactions occur when two or more mutations in different genes combine to result in a phenotype that is different from the expected phenotype when these mutations are tested separately in different individuals. 

For example, let A and B be the only two genes present in *P. pastoris* that code for an enzyme responsible of a limitant step in a vital pathway. When gene A is non-functional or missing due to a mutation, gene B can replace gene A's function and vice versa. This process allows *P. pastoris* to grow when a vital gene is deleted. Consequently, the phenotype when A or B are deleted separatedly in different individuals is **non-lethal**.

But if we delete gene A and B in the same individual, we will find that the phenotype is **lethal** because *P. pastoris* will not be able to grow up, because there is no gene this time that can replace the function of the lost vital genes.

Knowing all the above we can realise that gene A and B are interacting because when deleted in the same individual the obtained phenotype (lethal) is different than the obtained phenotype when the deletions occur in separated individuals. 

To determine how lethal the supression of genes can be, we measure the size of the colony resulting from an individual containing the suppressions that we want to study. The relation between the real and expected size of the colony is what we call **fitness**. Using this data and our mathematical model we are able to determine if there is an interaction between genes and what type of interaction is it.

## Types of genetic interaction
We will consider two main types of interaction:
1. **Negative genetic interaction**: Occurs when a combination of mutations leads to a fitness defect that is more exacerbated than expected. The previous genetic interaction example is in this type.
     - Synthetic lethality: Occurs when two non-letal mutations generate a non-viable mutant when combined.
2. **Positive genetic interaction**: Occurs when a combination of mutations leads to a fitness greater than expected.
     - Genetic suppression: Occurs when the mutations in the fitness defect of a query mutant is alleviated by a mutation in a second gene. 

## Algorithm
Our algorithm consists of two different coupled parts: training and testing. In the training we show a part of our data set to the algorithm so it is able to learn from it. In the testing we check how good our algorithm does its predictions about a part of the dataset that has not been never seen by the algorithm before.

In order to use this two parts the script must be run using the *fold(...)* method before doing anything. This will split our data set in five different chunks in which we will use [k-fold crossed validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)) in order to cross-validate the results from the algorithm.

After the split is done, you can run the algorithm using as input one of the data chunks that the method *fold(...)* has created to train the algorithm and its corresponding chunk to do the testing, consecutively.

After finishing the training, the algorithm will compute the testing results. Then, will store training and testing data from every sample in a separate csv file.

##### Training
The training part has many steps:
1. The algorithm gets the input data and digests it:
   - Every gene is a node in the network. 
   - Links are every assay between three genes and are tagged with 0 or 1 if there is interaction or not.
2. Many data structures are randomly initialized:
   - 2D-Matrix where every gene has its vector of possibilities of behaving like one of the groups of genes.
   - 3D-Matrix where every group of genes has a matrix of possibilities of interact with the other two group of genes.
3. The algorithm starts to iterate applying our model to maximize the **likelihood** parameter. This parameter describes how good our data model fits our experimental data.
4. When likelihood is high enough, the program will stop iterating and will save all the data from the model. 

The program will repeat step 2, 3 and 4 until the number of desired samples is completed.

##### Testing
After every sample has converged (its likelihood growing reached a certain threshold) testing results are computed for that sample. Testing and training results are stored in the same file for each sample.

## Usage
Our main script is "TrigenicInteractionPredictor.py". This file contains the class Model, that contains all the methods related to reading, processing and outputting data and also the methods related to training and testing. This script can only run one iteration at a time since there are data dependencies between interations. The script will iterate in each sample until likelihood converges. Precision, number of samples, number of iterations, likelihood frequency calculation, data source for training and other modifications to the algorithm behaviour can be applied using the defined parameters:

##### Informational Parameters

	-h, --help
		Displays the usage information.

##### Behaviour Parameters

	-i NUMITERATIONS, --iterations=NUMINTERATIONS
		Sets the maximum number of iterations that the program can compute for each sample. This does not mean that each sample is going to take NUMITERATIONS to finish, since the algorithm iterate over a sample until the likelihood reaches certain threshold. To force this behaviour, set frequencyCheck to 0, so likelihood is computed only at the end of NUMITERATIONS iterations and the program will not finish due to likelihood converge.

	-n NUMSAMPLES, --numSamples=NUMSAMPLES
		Sets how many samples will do the algorithm.

	-f LIKELIHOODFREQUENCYCHECK, --fcheck=LIKELIHOODFREQUENCYCHECK
		Sets how often the likelihood is computed. Every LIKELIHOODFREQUENCYCHECK iterations the likelihood will be computed. If LIKELIHOODFREQUENCYCHECK is set to 0, likelihood will be only computed at iteration number NUMITERATION, actually the last of the sample. If LIKELIHOODFREQUENCYCHECK is set to 1 likelihood will be computed every iteration. This behaviour parameter can be used with -b. LIKELIHOODFREQUENCYCHECK should be lower than NUMITERATIONS. Likelihood computing is computative and time expensive, so this parameter should be adjusted.

	-b LIKELIHOODCOMPUTINGSTARTINGITERATION, --bcheck=LIKELIHOODCOMPUTINGSTARTINGITERATION
		Sets in which iteration the algorithm will be able to start computing the likelihood. Can be used in conjunction with -f.

	-o OUTPUTPATH, --out=OUTPUTPATH
		Sets the output folder in which the program will store its result for every sample
	
##### Algorithm parameters

	-t TRAINFILEPATH, --train=TRAINFILEPATH
		Sets the path to the file containing the training information. This file should be generated using the *fold(...)* method.
	
	-e TESTFILEPATH, --test=TESTFILEPATH
		Sets the path to the file containing the testing information. This file should be generated using the *fold(...)* method.

	-k KLEVEL, --k=KLEVEL
		Sets the K argument to the algorithm. More K implies more computation but not necessarily better results or precision. Greater values of K may produce over-fitting. K should be 2 or greater tu use MMSBM. If not, algorithm is Single-Membership Stochastic Block Model.

##### Parallelization Parameters

	-s SAMPLENUMBERINITIALIZATION, --sampleIni=SAMPLENUMBERINITIALIZATION
		Sets in which number of sample we are currently. Used when running in parallel to be able to loop unroll the main loop of TrigenicInteractionPredictor.py and to avoid that different instances of the algorithm write in the same file. 
	
##### Default arguments
    --iterations=10000
    numSamples = 100
    sampleIni = 0
    frequencyCheck = 25
    beginCheck = 100
    train = '/home/aleixmt/Escritorio/TrigenicInteractionPredictor/data/folds/train0.dat'
    test = '/home/aleixmt/Escritorio/TrigenicInteractionPredictor/data/folds/test0.dat'
    outfilepath = ""
    argk = 2

## Code Health
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/51cacbf196634b1f81521e09bfdc9617)](https://www.codacy.com/app/AleixMT/TrigenicInteractionPredictor?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=AleixMT/TrigenicInteractionPredictor&amp;utm_campaign=Badge_Grade)

## Authors

* **Aleix Mariné** - [AleixMT](https://github.com/AleixMT) [aleix.marine@estudiants.urv.cat](aleix.marine@estudiants.urv.cat)
* **Marta Sales-Pardo** - [seeslab](https://github.com/seeslab) [marta.sales@urv.cat](marta.sales@urv.cat)
* **Roger Guimerà** - [roger.guimera@urv.cat](roger.guimera@urv.cat)
