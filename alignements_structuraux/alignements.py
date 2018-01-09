# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import numpy as np
from Bio.Seq import Seq
import Bio.Alphabet

seq1 = Seq("AGCATTTGGCTGGAAGCG")
print(seq1)
print(seq1.alphabet)
print(len(seq1))
seq2 = Seq("AGATGACTACCCTGGGTT")

def longer_common_sequence(seq1, seq2):
    m = np.zeros((len(seq1), len(seq2)))
    m[0, 0] = int(seq1[0] == seq2[0])
    for i in range(0, len(seq1)-1):
        for j in range(0, len(seq2)-1):
            m[i+1, j+1] = max(m[i+1, j], m[i, j+1], m[i, j] + int(seq1[i+1] == seq2[j+1]))
    return(m)
    
print(longer_common_sequence(seq1, seq2))