# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:19:05 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo

seq = Seq("WWAGCATTT-GGCTGG")

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
        self.seqs = self.getSeqs() + bloc.getSeqs()  #aligner(self.getSeqs, bloc)
        self.nbSeqs = self.getNbSeqs() + bloc.getNbSeqs()
        return 0 #retourner le score entre les deux
    
    def score(self, bloc):
        return 0
        
monBloc = bloc()

monBloc.init(seq)
print(monBloc.getNbSeqs())
print(monBloc.getSeqs())


score = monBloc.add(monBloc)
print(monBloc.getNbSeqs())
print(monBloc.getSeqs())
print(score)