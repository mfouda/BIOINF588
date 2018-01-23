# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:19:05 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo

seq = Seq("WWAGCATTTGGCTGG")

class bloc:
    
    def __init__(self):
        self.nbSeqs = 0
        self.seqs = []
        
    def init(self, seq):
        self.nbSeqs = 1
        self.seqs = [seq]
        
    def getNbSeqs(self):
        return self.nbSeqs
    
    def getSeqs(self):
        return self.seqs
    
    def getSeq(self, i):
        return self.getSeqs()[i]
    
    def add(self, bloc):
        self.seqs = aligner(self.getSeqs, bloc)
        self.nbSeqs = self.getNbSeqs() + bloc.getNbSeqs()
        
monBloc = bloc()

monBloc.init(seq)

print(monBloc.getNbSeqs())
print(monBloc.getSeqs())