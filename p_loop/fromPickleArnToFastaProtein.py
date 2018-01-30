# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:43:39 2018

@author: romai
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
import pickle as pk

dicoARN = pk.load(open( "dicoARN.p", "rb" ))
dicoPROT = dict()
for k, v in dicoARN.items():
    dicoPROT[k] = Seq(v, generic_rna).translate()
    
print(dicoARN["R0040C"])
print(dicoPROT["R0040C"])

pk.dump(dicoPROT, open( "dicoPROT.p", "wb" ) )


