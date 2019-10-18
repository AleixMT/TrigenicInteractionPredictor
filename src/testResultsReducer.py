#!/usr/bin/python3
# -*- coding: utf-8 -*-
#################################################################################################################
# Authors: Aleix Marine i Tena (aleix.marine@estudiants.urv.cat)                                                #
#       																									    #
# Description: Reads results from the results folder and computes the median, mean and standard deviation of    #
# the probability of interaction of every gene triplet separatedly for each fold and each K from all the test   #
# set.
# With this input_data, obtain plots of theresults
#################################################################################################################

import matplotlib
import getopt
import sys
import os
import codecs
import re
import copy
import math

if __name__ == "__main__":
    # Default arguments
    results_folder = "/home/aleixmt/Escritorio/TrigenicInteractionPredictor/data/DEFINITIVE_RESULTS/"

    # This block returns a value consisting of two elements: the first is a list of (option, value) pairs.
    # The second is the list of program arguments left after the option list was stripped.
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:", ["folder="])

        for opt, arg in opts:
            if opt in ("-f", "--folder"):
                if os.path.exists(str(arg)):
                    results_folder = arg
                else:
                    print("\n\nERROR: The selected path does not exist.")
                    raise ValueError

    except getopt.GetoptError:
        sys.exit(2)
    except ValueError:
        sys.exit(2)

    # input data structure
    input_data = {}

    # Data structure for each input triplete. triplete is array[2][4][5][].
    triplete_data = []
    for i in range(2):
        a = []
        for j in range(4):  # List of 4 positions corresponding to each K //HC
            b = []
            for k in range(5):  # List of 5 positions corresponding to each fold.
                # Notice that each triplete can't appear in all folds
                c = []  # Empty list to append probability in every sample.
                b.append(c)
            a.append(b)
        triplete_data.append(a)
    triplete_data[1] = []

    # Data structure for each output triplete
    triplete_result = []
    for i in range(4):   # List of 4 positions corresponding to each K //HC
        a = []
        for j in range(5):  # List of 5 positions corresponding to each fold.
            b = []  # Will be a list of 3 positions, corresponding to mean median and std. desviation
            a.append(b)
        triplete_result.append(a)
    likelihoods = copy.deepcopy(triplete_result)

    # find all the files recursively in the base folder
    print("· Start reading files ·")
    listGeneCodetoGeneName = []
    for i in range(5):
        listGeneCodetoGeneName.append([])

    print(listGeneCodetoGeneName)
    testCounter = 0
    for (dirpath, dirnames, filenames) in os.walk(results_folder):
        for f in filenames:
            file_pointer = os.path.join(dirpath, f)
            # Avoid and show possible casualties
            if os.stat(file_pointer).st_size == 0:  # Check if the file is empty
                print("· WARNING! Found empty file ·")
                continue
            else:
                if f[-1] == "#":  # check if it's a locked file
                    print("· WARNING! Found lock file. Skipping... ·")
                    continue  # and skip it if it is in

            # Get info of each sample using the path to the file (fold, sample number and K)
            k_number, fold_number, sample_number = file_pointer.split('/')[-3:]
            k_number = k_number.lstrip("K")
            fold_number = fold_number.lstrip("fold")
            sample_number = sample_number.split('_')[1]
            print("Dataset K=" + k_number + " fold=" + fold_number + " sample=" + sample_number)
            k_number = int(k_number) - 2  # Get the actual pointer in the list of Ks
            fold_number = int(fold_number)  # Explicit casting to int

            # Start reading the file
            # Begin with restoring a list that relates the gene number with the name of the gene
            file_handler = codecs.open(file_pointer)
            line = file_handler.readline()
            while not re.match("LIST OF REGISTERED GENES", line):
                line = file_handler.readline()
            file_handler.readline()  # Skip the line corresponding to the name of each column
            line = ""
            if not listGeneCodetoGeneName[fold_number]:  # We restore the list only one time per each fold
                line = file_handler.readline()  # Load first set of values
                line = line.rstrip('\n')  # strip carry return trailing character at end of the line
                while line:
                    fields = line.split("\t")
                    listGeneCodetoGeneName[fold_number].append(fields[1])
                    line = file_handler.readline()         # Load next set of values
                    line = line.rstrip('\n')  # strip carry return trailing character at end of the line

            # Read test result data and translate gene names
            file_handler = codecs.open(file_pointer)  # Reload file
            line = file_handler.readline()
            while not re.match("Held-out Likelihood", line):
                line = file_handler.readline()
            likelihoods[k_number][fold_number].append(float(line.split('\t')[1]))
            while not re.match("Test set:", line):  # Skip lines until find the string pattern "Test set:"
                line = file_handler.readline()
            file_handler.readline()  # Skip the line corresponding to the name of each column
            line = ""
            # Start reading input_data test
            while not line:
                line = file_handler.readline()
                print(line)
                line = line.rstrip('\n')  # strip carry return trailing character at end of the line
                if not line:  # If we find and empty string, corresponding to the end of the input_data block
                    break  # end reading
                interaction_probability, gene_triplet, real_interaction = line.split('\t')  # Get the data of each line
                interaction_probability = float(interaction_probability)
                real_interaction = int(real_interaction)

                if "307_532_78".__eq__(gene_triplet):
                    testCounter += 1
                    print(testCounter)
                geneNames = []
                for gene in gene_triplet.split("_"):
                    geneNames.append(listGeneCodetoGeneName[fold_number][int(gene)])

                geneNames.sort()
                geneNames = "_".join(geneNames)

                if geneNames in input_data.keys():
                    input_data[geneNames][0][k_number][fold_number].append(
                        interaction_probability)  # append interaction to the value in the dict
                else:  # if it's the first time you see this triplete
                    input_data[geneNames] = copy.deepcopy(triplete_data)  # CLONE data structure
                    input_data[geneNames][0][k_number][fold_number].append(
                        interaction_probability)  # append interaction to the value in the dict
                    input_data[geneNames][1].append(real_interaction)  # Append just one time

    print("· Reducing results ·")
    output = {}  # Output dictionary that relates triplete with its results (a copy of the variable triplete_results)
    # Calculate median, mean and std. dev. of probabilities of every triplet, for every fold in every k
    for key, value in input_data.items():
        sum_total = 0
        output[key] = copy.deepcopy(triplete_result)
        for k_number in range(len(value[0])):
            for fold_number in range(len(value[0][k_number])):
                # Every triplete appears only in one fold
                if value[0][k_number][fold_number]:  # enter if there are samples
                    print(key)
                    print(k_number)
                    print(fold_number)

                    print(value[0][k_number][fold_number])
                    # Mean
                    sum_total = 0
                    sum_real = 0
                    number_of_samples = len(value[0][k_number][fold_number])
                    print(number_of_samples)
                    for sample_number in range(number_of_samples):
                        sum_total += value[0][k_number][fold_number][sample_number]
                    mean = sum_total / number_of_samples
                    output[key][k_number][fold_number].append(mean)

                    # Median
                    value[0][k_number][fold_number].sort()
                    if number_of_samples % 2:  # if uneven
                        median = value[0][k_number][fold_number][round(number_of_samples / 2)]  # truncate
                    # Get the mean of the two elements in the center of the list
                    else:
                        half = int(number_of_samples / 2)
                        median = sum(value[0][k_number][fold_number][half - 1:half + 1]) / 2
                    output[key][k_number][fold_number].append(median)

                    # Std. Dev.
                    sum_square_diff = 0
                    for sample_number in range(number_of_samples):
                        sum_square_diff += (value[0][k_number][fold_number][sample_number] - mean) ** 2
                    stdDev = math.sqrt(sum_square_diff / len(value[0][k_number][fold_number]))
                    output[key][k_number][fold_number].append(stdDev)

    for k_number in range(len(likelihoods)):
        for fold_number in range(len(likelihoods[k_number])):
            mean_likelihood = float(sum(likelihoods[k_number][fold_number]) / len(likelihoods[k_number][fold_number]))
            likelihoods[k_number][int(fold_number)] = mean_likelihood

    print("· Remapping output ·")
    triplete_result_sorted = copy.deepcopy(triplete_result)
    print("ouputteee " + str(output))
    for key, value in output.items():
        for k_number in range(len(value)):
            for fold_number in range(len(value[k_number])):
                if value[k_number][fold_number]:

                    triplete_result_sorted[k_number][fold_number].append([key, value[k_number][fold_number][0], value[k_number][fold_number][1], value[k_number][fold_number][2], input_data[key][1][0]])

    print("· Compute metrics ·")
    positive_training_density = [0., 0., 0., 0., 0.]
    for fold_number in range(len(positive_training_density)):
        numPositives = 0
        counter = 0
        with open("../data/DATA_FOLDS/train" + str(fold_number) + ".dat", 'r') as f:
            for line in f.readlines():
                counter += 1
                if int(line.split('\t')[1]):
                    numPositives += 1
        positive_training_density[fold_number] = float(numPositives) / float(counter) # Obtain fraction of positives present in each training fold

    print(positive_training_density)
    print(triplete_result_sorted)
    for k_number in range(len(triplete_result_sorted)):
        for fold_number in range(len(triplete_result_sorted[k_number])):
            # AUC metric
            triplete_result_sorted[k_number][fold_number].sort(key=lambda tup: tup[2])  # Sort using the value of mean probability of interaction
            print(triplete_result_sorted[k_number][fold_number])
            predictedNumberOfPositivesTest = int(positive_training_density[fold_number] * len(triplete_result_sorted[k_number][fold_number]))
            counter = 0
            cutValue = 0

            print(triplete_result_sorted[0][0])
            for record in triplete_result_sorted[k_number][fold_number]:
                if predictedNumberOfPositivesTest == counter:
                    cutValue = record[1]
                    break
                counter += 1

            positives = []
            negatives = []
            counter = 0
            for record in triplete_result_sorted[k_number][fold_number]:
                if record[4]:  # record[4] is real interaction of the triplete
                    positives.append(record)
                else:
                    negatives.append(record)

            for positive in positives:
                for negative in negatives:
                    if positive[1] > negative[1]:  # Compare predicted probabilities
                        counter += 1
            auc = counter / (len(positives) * len(negatives))

            # Held-out likelihood
            heldoutLikelihood = 0.
            for sample in triplete_result_sorted[k_number][fold_number]:
                heldoutLikelihood += math.log(sample[1])  # likeli = sum(log(predicted)

            # precision, fallout, recall
            true_positives, false_positives, false_negatives, true_negatives = 0, 0, 0, 0
            for sample in triplete_result_sorted[k_number][fold_number]:
                predicted = sample[1]
                real = sample[4]
                if predicted >= cutValue:
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
            with open("/home/aleixmt/Escritorio/TrigenicInteractionPredictor/data/REDUCED_TRAIN_TEST_RESULTS/" + "K" + str(k_number + 2) + "_fold" + str(fold_number) + ".csv", 'w+') as f:
                f.write("\nHeld-OutLikelihood\tAUC\tPrecision\tRecall\tFallout\n")
                f.write(str(likelihoods[k_number][fold_number]) + '\t' + str(auc) + '\t' + str(precision) + '\t' + str(recall) + '\t' + str(fallout) + '\n')
                f.write("\nTripleteName\tMean\tMedian\tStdDev\tRealinteraction\n")
                for sample in triplete_result_sorted[k_number][fold_number]:
                    f.write(str(sample[0]) + '\t' + str(sample[1]) + '\t' + str(sample[2]) + '\t' + str(sample[3]) + '\t' + str(sample[4]) + "\n")

