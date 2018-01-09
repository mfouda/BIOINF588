# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import numpy as np
from Bio.Seq import Seq
import Bio.Alphabet

seq1 = Seq("AGCATTTGGCTGGAAGCG")

#print(seq1)
#print(seq1.alphabet)
#print(len(seq1))

seq2 = Seq("AGATGACTACCCTGGGTT")

def longer_common_sequence(seq1, seq2):
    m = np.zeros((len(seq1), len(seq2)))
    m[0, 0] = int(seq1[0] == seq2[0])
    
    for i in range(0, len(seq1)-1):
        m[i+1, 0] = max(m[i, 0], m[i, 0] + int(seq1[i+1] == seq2[0]))
        
    for j in range(0, len(seq2)-1):
        m[0, j+1] = max(m[0, j], m[0, j] + int(seq1[0] == seq2[j+1]))
            
    for i in range(0, len(seq1)-1):
        for j in range(0, len(seq2)-1):
            m[i+1, j+1] = max(m[i+1, j], m[i, j+1], m[i, j] + int(seq1[i+1] == seq2[j+1]))
            
    return(int(m[len(seq1) - 1, len(seq2) - 1]))
    
#print(longer_common_sequence(seq1, seq2))

def longer_common_sequence_arb(Seqs):
    print(len(Seqs))
    dims = []
    
    for i in range(0, len(Seqs)):
        dims += [len(Seqs[i])]
        
    m = np.zeros(dims)
    
    acc = np.zeros(len(Seqs))
    print(m)
    for s in range(0, len(Seqs)):
        acc = np.zeros(len(Seqs))
        #for i in range(0, len(Seqs[s])-1):
            #m[i+1, 0] = max(m[i, 0], m[i, 0] + int(seq1[i+1] == seq2[0]))
    
    return(m)
    
print(longer_common_sequence_arb([seq1, seq2]))
    