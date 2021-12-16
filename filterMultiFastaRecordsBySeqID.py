#!/common/opt/anaconda3/4.2.0/bin/python

import Bio
from Bio import SeqIO
import sys

##############################################################################
#USAGE filter.py <header.txt> <MultiFasta.fasta> <Outfile.fasta>             #
#Given a list of headers for a multifasta, this will pull those sequences    #
#out of the multifasta and put them in a specified outfile.                  #
##############################################################################

headerFile = sys.argv[1]
print(headerFile)
inFile = sys.argv[2]
print(inFile)
outFile = sys.argv[3]
print(outFile)

with open(headerFile) as f:
        headers = [x.strip('\n') for x in f]                  # strip line breaks from each line (x)

transcripts=[]                                                # initialize transcripts array
for i in headers:                                             # loop through each line
#    print(i)                                                 # uncomment to print each header for debugging
    transcripts.append(i.strip(">"))                          # add each line to the transcripts array, stripping the leading > if it's there

print("I have "+str(len(transcripts))+" records to find....") # Friendly chat script
print("Loading fasta to be filtered...")             

                                                              # This is the loop that does all the work
my_records=[]                                                 # initialize array to hold Seqio records
my_found_ids=[]                                               # initialize array to hold found ids
for seq_record in SeqIO.parse(inFile, "fasta"):               # parse the in file using SeqIO (each sequence becomes a record)
    #print(seq_record.id)                                     # uncomment for debugging to print each sequence id
    if seq_record.id in str(transcripts):                     # if the sequence id is in the array of transcript headers, then
        if seq_record.id not in my_found_ids:                 # if we haven't already found this id before (if it's not in found_ids)
            my_records.append(seq_record)                     # add it to the records array and 
            my_found_ids.append(seq_record.id)                # add it to the list of found ids
        else:                                                 # if we have already found it before, then 
            pass                                              # do nothing and go to the next record
#            print("Duplicate ? "+seq_record.id+" already found.") # uncomment for debugging
    else:                                                     # if it wasn't in the list of things to find, then 
#        #print(seq_record.id+" is not found :( ")            # uncomment for debugging
        pass                                                  # do nothing and move on 


for i in headers:                                             # finally, check to see if we found all the headers we were looking for 
    if i not in my_found_ids:                                 # if there's a header we were looking for that isn't in my_found_ids
        print(i+" not found?")                                # print that information to the screen 

print("I have found "+str(len(my_records))+ " records.")      # summarize in a friendly way 
SeqIO.write(my_records, outFile, "fasta")                     # This writes everything it found to the new fasta. 

