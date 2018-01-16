# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo

seq1 = Seq("AGCATTTGGCTGGAAGCG")

#print(seq1)
#print(seq1.alphabet)
#print(len(seq1))
print(Bio.SubsMat.MatrixInfo.blosum62['R', 'A'])
seq2 = Seq("AGATGACTACCCTGGGTT")

#seq1 = Seq("AAAAAAAAAAAAAAAAAA")
#seq2 = Seq("AAAAR")

def align(seq1, seq2, d):
    
    blosum62 = Bio.SubsMat.MatrixInfo.blosum62
    
    m = np.zeros((len(seq1) + 1, len(seq2) + 1)).astype(int)
    isfrom = np.zeros((len(seq1) + 1, len(seq2) + 1)).astype(int)
    
    for i in range(0, len(seq1)):     #Calcul de la première ligne
        m[i+1, 0] = m[i, 0] - d
        isfrom[i+1, 0] = 1
        
    for j in range(0, len(seq2)):     #Calcul de la première colonne
        m[0, j+1] = m[0, j] - d
        isfrom[0, j+1] = -1
    
    index, maxi = (0,0), m[0][0]
    for i in range(0, len(seq1)):     #Relation de récurrence
        for j in range(0, len(seq2)):
            
            if((seq1[i], seq2[j]) in blosum62):
                score = blosum62[seq1[i], seq2[j]]
            else:
                score = blosum62[seq2[j], seq1[i]]
            
            score += m[i, j]
            isfromij = 0
            
            if(m[i, j+1] - d > score):
                score = m[i, j+1] - d
                isfromij = 1
                
            if(m[i+1, j] - d > score):
                score = m[i, j+1] - d
                isfromij = -1
            
            isfrom[i+1, j+1] = isfromij
            m[i+1, j+1] = score
            
            if(score > maxi):
                maxi = score
                index = [i+1, j+1]
    
    print(m)
    print(isfrom)
    print(index, maxi)
    print(' ')
    
    str1 = ''
    str2 = ''
    
    while(index != (0, 0)):
        
        if(isfrom.item(tuple(index)) == 0):
            str1 = seq1[index[0] - 1] + str1
            str2 = seq2[index[1] - 1] + str2
            index = (index[0] - 1, index[1] - 1)
        
        if(isfrom.item(tuple(index)) == 1):
            str1 = seq1[index[0] - 1] + str1
            str2 = '-' + str2
            index = (index[0] - 1, index[1])
            
        if(isfrom.item(tuple(index)) == -1):
            str1 = '-' + str1
            str2 = seq2[index[1] - 1] + str2
            index = (index[0], index[1] - 1)
    
    print(seq1)
    print(seq2)
    print(' ')
    print(str1)
    print(str2)
    return('done')
    
align(seq1, seq2, 8)





    
# =============================================================================
#     
# def are_equals(Seqs, loc):      #Dit si tous les acides aminés sont les mêmes sur une location
#     b = True
#     for i in range(0, len(Seqs)):
#         for j in range(i + 1, len(Seqs)):
#             b = b and Seqs[i][loc[i]] == Seqs[j][loc[j]]
#     return b
# 
# def longest_common_sequence_arb(Seqs):
#     dims = []
#     
#     for i in range(0, len(Seqs)):
#         dims += [len(Seqs[i])]
#         
#     m = np.zeros(dims).astype(int)
#     
#     acc = tuple(np.zeros(len(Seqs)).astype(int))
#     m.itemset(acc, int(are_equals(Seqs, acc)))
#     
#     for s in range(0, len(Seqs)):
#         acc = np.zeros(len(Seqs)).astype(int)
#         for i in range(0, len(Seqs[s])-1):
#             maxi = m.item(tuple(acc))
#             acc[s] += 1
#             maxi = max(maxi, maxi + int(are_equals(Seqs, acc)))
#             m.itemset(tuple(acc), maxi)
#             
#     boucles = []
#     
#     for s in range(0, len(Seqs)):
#         boucles += [[1, len(Seqs[s])]]
#     
#     # initialisation de la boucle
#     boucles = [(type(x)==list and [x] or [[x]])[0] for x in boucles]
#     boucles = [(len(x)==1 and [[0] + x] or [x])[0] for x in boucles]
#     boucles = [(len(x)==2 and [x + [1]] or [x])[0] for x in boucles]
#     compteurs = [x[0] for x in boucles]  # => initialisation des compteurs de boucle (= valeur initiale de boucle)
#     compteurs[-1] -= boucles[-1][2]  # le 1er incrément sera pour rien
#     indmax = len(boucles)-1  # => index de la dernière boucle (la plus interne)
#     finboucle = False  # => drapeau pour condition de fin de la boucle globale
#      
#     while True:
#         # incrémentation des compteurs, test de boucle et condition de sortie
#         for x in range(indmax, -1, -1):
#             compteurs[x] += boucles[x][2]
#             if compteurs[x]>=boucles[x][1]:
#                 if x==0:
#                     finboucle = True
#                     break
#                 compteurs[x] = boucles[x][0]
#             else: break
#         if finboucle: break
#         
#         subcount = compteurs.copy()
#         subsubcount = compteurs.copy()
#         subcount[0] -= 1
#         subsubcount[0] -= 1
#         maxi = m.item(tuple(subcount))
#         
#         for i in range(1, len(Seqs)):
#             subcount = compteurs.copy()
#             subcount[i] -= 1
#             subsubcount[i] -= 1
#             maxi = max(maxi, m.item(tuple(subcount)))
#         
#         maxi = max(maxi, m.item(tuple(subsubcount)) + int(are_equals(Seqs, compteurs)))
#         m.itemset(tuple(compteurs), maxi)
#     
#     return(m.item(tuple([-1]*len(Seqs))))
#     
# #print("longest_common_sequence_arb =",longest_common_sequence_arb([seq1, seq2, seq2]))
#     
# =============================================================================
