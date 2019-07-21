#!/usr/bin/python3
# -*- coding: utf-8 -*-
#################################################################################################################
# Authors: Aleix Marine i Tena (aleix.marine@estudiants.urv.cat)                                                #
#       																									    #
# Description: Reads results data from the results folder and computes the median, mean and standard deviation  #
# of the probabilities of interaction for every gene triplet in the test set.                                   #
# With this data, obtain plots of theresults
#################################################################################################################

import matplotlib.pyplot as plt
import getopt
import sys
import os
import codecs
import re

if __name__ == "__main__":
    # Default arguments
    results_folder = "/home/aleixmt/DEFINITIVE_RESULTS"

    # This method returns value consisting of two elements: the first is a list of (option, value) pairs.
    # The second is the list of program arguments left after the option list was stripped.
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:", ["folder="])

        for opt, arg in opts:
            if opt in ("-f", "--folder"):  # Show the usage if help is called
                if os.path.exists(str(arg)):
                    results_folder = arg
                else:
                    print("\n\nERROR: The selected path does not exist.")
                    raise ValueError

    except getopt.GetoptError:
        sys.exit(2)
    except ValueError:
        sys.exit(2)
    data = {}
    for (dirpath, dirnames, filenames) in os.walk(results_folder):
        for f in filenames:
            file_pointer = os.path.join(dirpath, f)
            if os.stat(file_pointer).st_size == 0:  # Check if the file is empty
                continue
            else:
                pass
            k_number, fold_number, sample_number = file_pointer.split('/')[-3:]  # Get the data form the last part of the path
            sample_number = sample_number.split('_')[1]
            print("Dades fitxer: " + k_number + fold_number + sample_number)
            file_handler = codecs.open(file_pointer)
            line = file_handler.readline()

            while not re.match("Test set:", line):  # Skip lines until find the string pattern "Test set:"
                line = file_handler.readline()
            file_handler.readline()  # Skip the line corresponding to the name of each column
            for line in file_handler.readlines():  # Start reading data test
                line = line.rstrip()  # strip carry return trailing character at end of the line
                if not line:  # If we find and empty string, corresponding to the end of the data block
                    break  # end reading
                interaction_probability, gene_triplet, real_interaction = line.split('\t')
                interaction_probability = float(interaction_probability)
                try:
                    data[gene_triplet].append(interaction_probability)  # append interaction to the value in the dict
                except KeyError:  # if it's the first time you see this triplete
                    data[gene_triplet] = [interaction_probability]

    result_data = {}
    for key, value in data.items():
        sum_total = 0
        for i in value:
            sum_total += i
        mean = sum_total / len(value)
        result_data[key] = [mean]
        if len(value) % 2:
            median = value[round(len(value) / 2)]
        else:
            median = (value[round(len(value) / 2)] + value[round(len(value) / 2) + 1]) / 2
        result_data[key].append(median)

    print("\n\n" + str(result_data))










# Contains all methods and data structures of the algorithm.