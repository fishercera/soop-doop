#!/common/opt/anaconda3/4.2.0/bin/python

import Bio
from Bio import SeqIO
import sys

##############################################################################
#USAGE filter.py <header.txt> <MultiFasta.fasta> <Outfile.fasta>             #
#Given a list of headers for a multifasta, this will skip those sequences    #
#from the specificed multifasta and put the rest into the outfile.           #
##############################################################################

headerFile = sys.argv[1]
print(headerFile)
inFile = sys.argv[2]
print(inFile)
outFile = sys.argv[3]
print(outFile)

with open(headerFile) as f:
	headers = [x.strip('\n') for x in f]

transcripts=[]
for i in headers:
#    print(i)
    transcripts.append(i.strip(">"))

print("I have "+str(len(transcripts))+" records to find....")
print("Loading fasta to be filtered...")

#transcript=["transcrip1", "transcript2", "transcript3"...."transcriptn"]
#TODO add sys.arg to read through a list of headers
#TODO add sys.arg to get master fasta file

my_records=[]
my_kept_ids=[]
for seq_record in SeqIO.parse(inFile, "fasta"):
    #print(seq_record.id)
    if seq_record.id not in str(transcripts): ### Only keeping the ones not found in the header file.
        if seq_record.id not in my_kept_ids:
            my_records.append(seq_record)
            my_kept_ids.append(seq_record.id)
        else:
            print("Duplicate ? "+seq_record.id+" already found.")
    else:
        #print(seq_record.id+" is not found :( ")
        pass


print("I have found "+str(len(my_records))+ " records.")
SeqIO.write(my_records, outFile, "fasta")
