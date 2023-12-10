#!/usr/bin/python

import Bio
from Bio import SeqIO
import sys

##############################################################################
#USAGE: given the ID of a sequence, fetches that sequence and writes it      #
#to a file with a gene-name appended.                                        #
# getSeq.py <assembly.fasta> <seqID> <gene-name>                             #
#OUTPUT: seqID.gene-name.fasta                                               #
##############################################################################

#headerFile = sys.argv[1]
#print(headerFile)

inFile = sys.argv[1]
print(inFile)
seqID = sys.argv[2]
print(seqID)
geneName = sys.argv[3]

header = seqID.strip(">")
print(header)

outFile = header+"."+geneName+".fasta"
print(outFile)

# This is the loop that does all the work
my_records=[]
assembly =  SeqIO.to_dict(SeqIO.parse(inFile, "fasta"))

for key in assembly.keys():
    if key.startswith(header):
        print(key)
        print(assembly[key])
        with open(outFile, "a") as output:
            SeqIO.write(assembly[key], output, "fasta")
    else:
        pass

print(len(assembly))
