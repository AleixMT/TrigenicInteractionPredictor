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
    results_folder = "/media/_Shared/DEFINITIVE_RESULTS"

    # This method returns value consisting of two elements: the first is a list of (option, value) pairs.
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

    # find all the files recursively in the base folder
    print("· Start reading files ·")
    for (dirpath, dirnames, filenames) in os.walk(results_folder):
        for f in filenames:
            file_pointer = os.path.join(dirpath, f)
            if os.stat(file_pointer).st_size == 0:  # Check if the file is empty
                print("Empty file")
                continue
            else:
                if f[-1] == "#":  # check if it's a locked file
                    print("Found lock file. Skipping...")
                    continue  # and skip it if it is in
                pass
            #Get info of each sample using the path to the file
            k_number, fold_number, sample_number = file_pointer.split('/')[-3:]  # Get the input_data form the last part of the path
            k_number = k_number.lstrip("K")
            fold_number = fold_number.lstrip("fold")
            sample_number = sample_number.split('_')[1]
            print("Dataset K=" + k_number + " fold=" + fold_number + " sample=" + sample_number)
            k_number = int(k_number) - 2  # Get the actual pointer in the list of Ks
            fold_number = int(fold_number)  # Explicit casting to int

            # Start reading the file
            file_handler = codecs.open(file_pointer)
            line = file_handler.readline()
            while not re.match("Test set:", line):  # Skip lines until find the string pattern "Test set:"
                line = file_handler.readline()
            file_handler.readline()  # Skip the line corresponding to the name of each column

            # Start reading input_data test
            for line in file_handler.readlines():
                #print(line)  # Testing
                line = line.rstrip('\n')  # strip carry return trailing character at end of the line
                if not line:  # If we find and empty string, corresponding to the end of the input_data block
                    break  # end reading
                interaction_probability, gene_triplet, real_interaction = line.split('\t')  # Get the data of each line
                #print("Gene triplete" + gene_triplet + " has " + interaction_probability + "prob. ")

                interaction_probability = float(interaction_probability)
                real_interaction = int(real_interaction)
                if gene_triplet in input_data:
                    input_data[gene_triplet][0][k_number][fold_number].append(interaction_probability)  # append interaction to the value in the dict
                    input_data[gene_triplet][1].append(real_interaction)  # TESTING
                else:  # if it's the first time you see this triplete
                    input_data[gene_triplet] = copy.deepcopy(triplete_data)  # CLONE data structure
                    input_data[gene_triplet][0][k_number][fold_number].append(interaction_probability)  # append interaction to the value in the dict
                    input_data[gene_triplet][1].append(real_interaction)  # TESTING


    print(input_data["305_804_907"][1])  # TESTING
    '''
    input()
    #print(input_data)
    #TESTING
    for key, value in input_data.items():
        for k_number in range(len(value[0])):
            for fold_number in range(5):
                # print("sum real number: " + str(value[1][k_number][fold_number]))
                output_filename = "K" + str(k_number + 2) + "_fold" + str(fold_number)
                if os.path.exists(output_filename):
                    with open(output_filename, 'a') as f:
                        f.write(str(key) + "\t")
                        for prob in value[0][k_number][fold_number]:
                            f.write(str(prob) + "\t")
                        f.write("\n")
                else:
                    with open(output_filename, 'a+') as f:
                        f.write("Gene triplete\tProb\n")
                        f.write(str(key) + "\t")
                        for prob in value[0][k_number][fold_number]:
                            f.write(str(prob) + "\t")
                        f.write("\n")
    # TESTING
    '''
    output = {}
    print("· Computing Results ·")
    # Calculate median, mean and std. dev. of probabilities of every triplet, for every fold in every k
    for key, value in input_data.items():
        sum_total = 0
        output[key] = copy.deepcopy(triplete_result)
        for k_number in range(len(value[0])):
            for fold_number in range(len(value[0][k_number])):
                # Every triplete appears only in one fold
                if value[0][k_number][fold_number]:  # enter if there are samples
                    # Mean
                    sum_total = 0
                    sum_real = 0
                    number_of_samples = len(value[0][k_number][fold_number])
                    for sample_number in range(number_of_samples):
                        sum_total += value[0][k_number][fold_number][sample_number]
                    mean = sum_total / number_of_samples
                    output[key][k_number][fold_number].append(mean)

                    # Median
                    value[0][k_number][fold_number].sort()
                    if number_of_samples % 2:  # if not even
                        median = value[0][k_number][fold_number][round(number_of_samples / 2)]
                    else:
                        half = int(number_of_samples / 2)
                        median = sum(value[0][k_number][fold_number][half - 1:half])
                        print(value[0][k_number][fold_number][half - 1:half])  # TESTING
                    output[key][k_number][fold_number].append(median)

                    # Std. Dev.
                    sum_square_diff = 0
                    for sample_number in range(number_of_samples):
                        sum_square_diff += (value[0][k_number][fold_number][sample_number] - mean) ** 2
                    stdDev = math.sqrt(sum_square_diff / (len(value[0][k_number][fold_number]) - 1))
                    output[key][k_number][fold_number].append(stdDev)

    print("· Writing output ·")
    os.system("rm K*")
    for key, value in output.items():
        for k_number in range(len(value[0])):
            for fold_number in range(5):
                #print("sum real number: " + str(value[1][k_number][fold_number]))
                if value[k_number][fold_number]:
                    output_filename = "K" + str(k_number) + "_fold" + str(fold_number)
                    if os.path.exists(output_filename):
                        with open(output_filename, 'a') as f:
                            # Key   Mean    Median  StdDev  Real
                            f.write(str(key) + "\t" + str(value[k_number][fold_number][0]) + "\t" + str(value[k_number][fold_number][1]) + "\t" + str(value[k_number][fold_number][2]) + "\t" + str(input_data[key][1].pop()) + "\n")
                    else:
                        with open(output_filename, 'a+') as f:
                            f.write("Triplete\tMean\tMedian\tStdDev\tRealinteraction\n")
                            f.write(str(key) + "\t" + str(value[k_number][fold_number][0]) + "\t" + str(value[k_number][fold_number][1]) + "\t" + str(value[k_number][fold_number][2]) + "\t" + str(input_data[key][1].pop()) + "\n")

