#!/usr/bin/python3
# -*- coding: utf-8 -*-
#################################################################################################################
# Authors: Aleix Marine i Tena (aleix.marine@estudiants.urv.cat)                                                #
#       																									    #
# Description: Reads results input_data from the results folder and computes the median, mean and standard deviation  #
# of the probabilities of interaction for every gene triplet in the test set.                                   #
# With this input_data, obtain plots of theresults
#################################################################################################################

import matplotlib
import getopt
import sys
import os
import codecs
import re

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

    # Data structures
    input_data = {}
    output = {}

    # Data structure for each input triplete. triplete is array[2][4][5][100]
    triplete_data = []
    for i in range(2):
        a = []
        for j in range(4):  # List of 4 positions corresponding to K //HC
            b = []
            for k in range(5):  # List of 5 positions corresponding to each fold
                c = []  # List of 100 positions corresponding to samples //HC
                b.append(c)
            a.append(b)
        triplete_data.append(a)

    # Data structure for each output triplete
    triplete_result = []
    for i in range(4):
        a = []
        for j in range(5):
            b = [0 for _ in range(2)]
            a.append(b)
        triplete_result.append(a)

    # find all the files recursively in the base folder
    for (dirpath, dirnames, filenames) in os.walk(results_folder):
        for f in filenames:
            file_pointer = os.path.join(dirpath, f)
            if os.stat(file_pointer).st_size == 0:  # Check if the file is empty
                print("Empty file")
                continue
            else:
                if f[-1] == "#":  # check if it's a locked file
                    print("Found lock file. Skipping...")
                    continue  # and skip it if it is
                pass
            #Get info of each sample using the path to the file
            k_number, fold_number, sample_number = file_pointer.split('/')[-3:]  # Get the input_data form the last part of the path
            k_number = k_number.lstrip("K")
            fold_number = fold_number.lstrip("fold")
            sample_number = sample_number.split('_')[1]
            print("Dataset K=" + k_number + " fold=" + fold_number + " sample=" + sample_number)

            # Start reading the file
            file_handler = codecs.open(file_pointer)
            line = file_handler.readline()
            while not re.match("Test set:", line):  # Skip lines until find the string pattern "Test set:"
                line = file_handler.readline()
            file_handler.readline()  # Skip the line corresponding to the name of each column

            # Start reading input_data test
            for line in file_handler.readlines():
                line = line.rstrip('\n')  # strip carry return trailing character at end of the line
                if not line:  # If we find and empty string, corresponding to the end of the input_data block
                    break  # end reading
                interaction_probability, gene_triplet, real_interaction = line.split('\t')
                interaction_probability = float(interaction_probability)
                real_interaction = int(real_interaction)

                try:
                    input_data[gene_triplet][0][int(k_number) - 2][int(fold_number)].append(interaction_probability)  # append interaction to the value in the dict

                except KeyError:  # if it's the first time you see this triplete
                    input_data[gene_triplet] = triplete_data.copy()  # Create data structure
                    input_data[gene_triplet][0][int(k_number) - 2][int(fold_number)].append(interaction_probability)  # append interaction to the value in the dict
                input_data[gene_triplet][1][int(k_number) - 2][int(fold_number)].append(real_interaction)  # append interaction to the value in the dict

    print("Data reading finished")
    # Calculate median and mean of probabilities of every triplet, for every fold in every k
    for key, value in input_data.items():
        sum_total = 0
        output[key] = triplete_result.copy()
        for k_number in range(len(value[0])):
            for fold_number in range(5):
                sum_total = 0
                sum_real = 0
                print("key " + str(key) + " k " + str(k_number) + " fold " + str(fold_number) + " has:" + str(len(value[0][k_number][fold_number])) + " samples")
                for sample_number in range(len(value[0][k_number][fold_number])):
                    sum_total += value[0][k_number][fold_number][sample_number]
                    sum_real += value[1][k_number][fold_number][sample_number]
                mean = sum_total / len(value[0][k_number][fold_number])
                output[key][k_number][fold_number][0] = mean
                print("i té de suma " + str(sum_real))
                #value[1][k_number][fold_number] = sum_real
                value[0][k_number][fold_number].sort()
                number_of_samples = len(value[0][k_number][fold_number])
                if number_of_samples % 2:
                    median = value[0][k_number][fold_number][round(number_of_samples / 2)]
                else:
                    median = sum(value[0][k_number][fold_number][round(number_of_samples / 2) - 1:round(number_of_samples / 2)]) / 2
                output[key][k_number][fold_number][1] = median

    print("writing output")
    for key, value in output.items():
        for k_number in range(len(value[0])):
            for fold_number in range(5):
                #print("sum real number: " + str(value[1][k_number][fold_number]))
                output_filename = "K" + str(k_number) + "_fold" + str(fold_number)
                if os.path.exists(output_filename):
                    with open(output_filename, 'a') as f:
                        f.write(str(value[k_number][fold_number][0]) + "\t" + str(key) + "\t" + str(input_data[key][1][k_number][fold_number]) + "\n")
                else:
                    with open(output_filename, 'a+') as f:
                        f.write("Mean Interaction Probability\tGene triplete\tReal interaction\n")
                        f.write(str(value[k_number][fold_number][0]) + "\t" + str(key) + "\t" + str(input_data[key][1][k_number][fold_number]) + "\n")

