# TrigenicInteractionPredictor
Algorithm that predicts interaction between triplets of genes.

## Description

Genetic interactions occur when mutations in different genes combine to result in a phenotype that is different from expected when observed in different individuals. In other words, when the phenotype of an individual that has two different mutations is significable different from the phenotype of an individual with one or the other mutation individiually, it has occured a genetic interaction.
* This algorithm applies the model Mixed-Membership [Stochastic Block Model](https://en.wikipedia.org/wiki/Stochastic_block_model) (MMSBM) to predict interaction between tripletes of genes in [Pichia pastoris](https://en.wikipedia.org/wiki/Pichia_pastoris). 
* We use supplementary materials from the article [Systematic analysis of complex genetic interactions"](http://science.sciencemag.org/content/360/6386/eaao1729). DOI: 10.1126/science.aao1729 to get our data.
* These materials can be found [here](https://www.dropbox.com/sh/4wblbdwzy4bki53/AACM46GqkfJmzS7iekKcG4Wba?dl=0). Nevertheless, they will move soon to git lfs in order to have the datasets more available.

## Algorithm

The main program has many steps:
1.- The algorithm gets the input data and digests it:
  * Every gene is a node in the network. 
  * Links are every assay between three genes.
    * This links are tagged with 0 or 1 if there is interaction or not.  
2.- Many data structures are initialized with the input:
  * 2D-Matrix where every gene has its vector of possibilities of behaving like one of the groups of genes.
  * 3D-Matrix where every group of genes has a matrix of possibilities of interact with the other two group of genes.
3.- The algorithm starts to iterate to get the maximization of the likelihood.
4.- If gets the covergement a message will be shown.
5.- All the data from the execution will be saved automatically for every sample even if it doesn't converge.

## Usage

To execute the code and speeding it up you'll use a **pypy3** virtual environment. There's one uploaded in the repository in the folder `src/`. To execute it open a terminal and situate it in the `src` folder. Then type, for executing the algorithm with the default options:
```
pypy3-v6.0.0-linux64/bin/pypy3 TrigenicInteractionPredictor.py
```
These default arguments are:
* iterations = 1000
* samples = 1
* frequencyCheck = 1
* filename = "Data_S1.csv"
* interactionType = "trigenic"
* cutOffValue = -0.08
* argk = 10

You can use other arguments, specifying **_all_** of them after the text that executes the program:
1.- iterations[positive integer]: Number of iterations done by algorithm.
2.- samples[positive integer]: Number of samples done by algorithm.
3.- frequencyCheck[positive integer] Number of iterations needed to check if likelihood has converged.
4.- filename[string]: Name of the dataset filename.
5.- interactionType[string:{Trigenic,Digenic,\*}]: Type of interaction selected.
6.- cutOffValue[real]: Value used to determine if an interaction is positive or negative.
7.- argk[integer]: Number of groups to use in the algorithm (Increases lineally the computation cost).

For example:
```
pypy3-v6.0.0-linux64/bin/pypy3 TrigenicInteractionPredictor.py 1000 5 10 Data_S1.csv Trigenic -0.08
```

## Code Health
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/51cacbf196634b1f81521e09bfdc9617)](https://www.codacy.com/app/AleixMT/TrigenicInteractionPredictor?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=AleixMT/TrigenicInteractionPredictor&amp;utm_campaign=Badge_Grade)

## Authors

* **Aleix Mariné** - [AleixMT](https://github.com/AleixMT) [aleix.marine@estudiants.urv.cat](aleix.marine@estudiants.urv.cat)
* **Marta Sales-Pardo** - [seeslab](https://github.com/seeslab) [marta.sales@urv.cat](marta.sales@urv.cat)
* **Roger Guimerà** - [roger.guimera@urv.cat](roger.guimera@urv.cat)
