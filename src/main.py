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
# first argument. Data from the first file is then written in UTF-8 format with the same name as the input file #
# plus a 'utf-8.csv' suffix. If output file exists, its content will be overwritten.                            #
#                                                                                                               #
# -Argument 1: relative route from the same folder where the script is executed.                                #
#################################################################################################################

import re
import codecs
import sys
    
from enum import Enum

class Gene:
    
    def __init__(self, name, interaction):
        self.name = name
        self.interactions = [interaction]
        self.vector = []
        
    def __init__(self, name, strain):
        self.name = name
        self.interactions = []
        self.vector = []
        self.strain = strain
        
    def addInteraction(self, interaction):
        self.interactions.append(interaction)
        
class Group:
    
    def __init__(self):
        genes = []
        
    def __init__(self, genePointer):
        genes = [genePointer]
        
class Pointer:
        
    def __init__(self, value):
        if type(self) is Pointer:
            raise Exception('Base is an abstract class and cannot be instantiated directly')
        
        self.value = value

class GenePointer(Pointer):
    
    def __init__(self, gene, value):
        super().__init__(value)
        self.gene = gene

class GroupPointer(Pointer):
    
    def __init__(self, group, value):
        super().__init__(value)
        self.group = group
        
class Interaction:
    
    def __init__(self, genes, score):
        if type(self) is Interaction:
            raise Exception('Base is an abstract class and cannot be instantiated directly')
        
        self.genes.extends(genes)
        self.score = score
        
class DigenicInteraction(Interaction):

    def __init__(self, genes, interactionType):
        self.interactionType = interactionType
        super().__init__(genes)
        
class TrigenicInteraction(Interaction):

    def __init__(self, genes, interactionType):
        self.interactionType = interactionType
        super().__init__(genes)
        
class TypeTrigenic(Enum):
    novel = auto()
    unclassified = auto()
    modifiedQm = auto()
    modifiedQmAm = auto()
    modifiedQmAM = auto()
    modifiedQmAMm = auto()
    modifiedAm = auto()
    modifiedAM = auto()
    modifiedAmM = auto()
    
class TypeDigenic(Enum):
    digenic = auto()
    
def getInput(genes, interaction):
    dataset = codecs.open(sys.argv[1], encoding='utf-8', mode='r+')
    for line in dataset:
        fields = re.split(r'\t+', line)
        querystrains = fields[0].split('+')
        queryalleles = fields[1].split('+')
        arraystrain = fields[2]
        arrayallele = fields[3]
        interactionType = fields[4]
        score = fields[5]
        Pvalue = fields[6]
        for gene 
        if queryalleles[0] in genes:
            genes.
        
    
    return 0
    
if __name__ == "__main__":
    
    getInput([], [])
    


    

