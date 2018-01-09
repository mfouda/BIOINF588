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
    m = np.zeros((len(seq1), len(seq2))).astype(int)
    m[0, 0] = int(seq1[0] == seq2[0])
    
    for i in range(0, len(seq1)-1):
        m[i+1, 0] = int(max(m[i, 0], m[i, 0] + int(seq1[i+1] == seq2[0])))
        
    for j in range(0, len(seq2)-1):
        m[0, j+1] = int(max(m[0, j], m[0, j] + int(seq1[0] == seq2[j+1])))
            
    for i in range(0, len(seq1)-1):
        for j in range(0, len(seq2)-1):
            m[i+1, j+1] = int(max(m[i+1, j], m[i, j+1], m[i, j] + int(seq1[i+1] == seq2[j+1])))
    print(m)
    return(int(m[len(seq1) - 1, len(seq2) - 1]))
    
print("longer_common_sequence =",longer_common_sequence(seq1, seq2))
    
    
def are_equals(Seqs, loc):
    b = True
    for i in range(0, len(Seqs)):
        for j in range(i + 1, len(Seqs)):
            b = b and Seqs[i][loc[i]] == Seqs[j][loc[j]]
    return b

def longer_common_sequence_arb(Seqs):
    dims = []
    
    for i in range(0, len(Seqs)):
        dims += [len(Seqs[i])]
        
    m = np.zeros(dims).astype(int)
    
    acc = tuple(np.zeros(len(Seqs)).astype(int))
    m.itemset(acc, int(are_equals(Seqs, acc)))
    
    for s in range(0, len(Seqs)):
        acc = np.zeros(len(Seqs)).astype(int)
        for i in range(0, len(Seqs[s])-1):
            maxi = m.item(tuple(acc))
            acc[s] += 1
            maxi = max(maxi, maxi + int(are_equals(Seqs, acc)))
            m.itemset(tuple(acc), maxi)
            
    boucles = []
    
    for s in range(0, len(Seqs)):
        boucles += [[1, len(Seqs[s])]]
    
    # initialisation de la boucle
    boucles = [(type(x)==list and [x] or [[x]])[0] for x in boucles]
    boucles = [(len(x)==1 and [[0] + x] or [x])[0] for x in boucles]
    boucles = [(len(x)==2 and [x + [1]] or [x])[0] for x in boucles]
    compteurs = [x[0] for x in boucles]  # => initialisation des compteurs de boucle (= valeur initiale de boucle)
    compteurs[-1] -= boucles[-1][2]  # le 1er incrément sera pour rien
    indmax = len(boucles)-1  # => index de la dernière boucle (la plus interne)
    finboucle = False  # => drapeau pour condition de fin de la boucle globale
     
    while True:
        # incrémentation des compteurs, test de boucle et condition de sortie
        for x in range(indmax, -1, -1):
            compteurs[x] += boucles[x][2]
            if compteurs[x]>=boucles[x][1]:
                if x==0:
                    finboucle = True
                    break
                compteurs[x] = boucles[x][0]
            else: break
        if finboucle: break
    
        
        #maxi = m.item(tuple(acc))
        #acc[s] += 1
        #maxi = max(maxi, maxi + int(are_equals(Seqs, acc)))
        #m.itemset(tuple(acc), maxi)
        
        subcount = compteurs.copy()
        subsubcount = compteurs.copy()
        subcount[0] -= 1
        subsubcount[0] -= 1
        maxi = m.item(tuple(subcount))
        for i in range(1, len(Seqs)):
            subcount = compteurs.copy()
            subcount[i] -= 1
            subsubcount[i] -= 1
            maxi = max(maxi, m.item(tuple(subcount)))
        
        maxi = max(maxi, m.item(tuple(subsubcount)) + int(are_equals(Seqs, compteurs)))
        m.itemset(tuple(compteurs), maxi)
    
    return(m.item(tuple([-1]*len(Seqs))))
    
print("longer_common_sequence_arb =",longer_common_sequence_arb([seq1, seq2, seq2]))
    