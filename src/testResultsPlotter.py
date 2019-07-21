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

    for (dirpath, dirnames, filenames) in os.walk(results_folder):
        for f in filenames:
            file_pointer = os.path.join(dirpath, f)
            k_number = file_pointer.split('/')[-3]
            fold_number = file_pointer.split('/')[-2]
            sample_number = file_pointer.split('/')[-1].split('_')[1]
            print(k_number + fold_number + sample_number)



# Contains all methods and data structures of the algorithm.