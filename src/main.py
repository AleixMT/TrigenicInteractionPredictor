#!/usr/bin/python
import re

def getInput():
    print "cosa"
    with open('aao1729_Data_S2.txt', 'r+') as input:
        for line in input:
            columns = re.split(r'\t+', line)
            print columns

    return 0

#if __name__ == "__main__":
getInput()
print "jelou"


    

    

  