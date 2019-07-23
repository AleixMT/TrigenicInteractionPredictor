#!/usr/bin/python

#################################################################################################################
# Authors: Aleix Marine i Tena (aleix.marine@estudiants.urv.cat)                                                #
# Last review: 27/6/2018                                                                                        #
# Version: 1.0                                                                                                  #
#                                                                                                               #
# Permissions: Needs permissions to read file in the same folder as the script, specified by the first argument,#
# and permissions to create a new file in the folder of the script.                                             #
#                                                                                                               #
# Parameters and descriptions: Reads the UTF-16 file in the folder where the script is located specified by the #
# first argument. Data from the first file is then written in UTF-8 format with the same name as the input_data file #
# plus a 'utf-8.csv' suffix. If output file exists, its content will be overwritten.                            #
#                                                                                                               #
# -Argument 1: relative route from the same folder where the script is executed.                                #
#################################################################################################################

import sys
import codecs

try:
    datasetutf16 = codecs.open(sys.argv[1], encoding='utf-16', mode='r+')
    datasetutf8 = codecs.open(sys.argv[1]+'utf8.csv', encoding='utf-8', mode='w+')
    datasetutf8.write(datasetutf16.read())
    datasetutf8.close()
    datasetutf16.close()
except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)