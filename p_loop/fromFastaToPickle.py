# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:43:39 2018

@author: romai
"""

from Bio import SeqIO
import pickle as pk

input_file = "data/orf_coding_all.fasta"
fasta_sequences = SeqIO.parse(open(input_file),'fasta')

dicoARN = dict()
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    dicoARN[name] = sequence
    
pk.dump(dicoARN, open( "dicoARN.p", "wb" ) )
#print(dicoARN)