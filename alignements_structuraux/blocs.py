# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:19:05 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo

class bloc:
    
    def __init__(self, seq):
        self.nbSeqs = 1
        self.seqs = [seq]
        self.score = np.nan
        
    def getNbSeqs(self):
        return self.nbSeqs
    
    def getSeqs(self):
        return self.seqs
    
    def getSeq(self, i):
        return self.getSeqs()[i]
    
    def getCol(self, i):
        col = []
        for k in range(0, self.getNbSeqs()):
            col += [self.getSeq(k)[i]]
        return col
    
    def show(self):
        print('############')
        print('### BLOC ###')
        print('############')
        if(self.getNbSeqs() == 1):
            print('The bloc has 1 sequence')
        else:
            print('The bloc has '+str(self.getNbSeqs())+' sequences')
        if(not np.isnan(self.score)):
            print('The last merging score is',self.score)
        print('-'*len(self.getSeq(0)))
        for i in range(0, self.getNbSeqs()):
            print(self.getSeq(i))
        print('-'*len(self.getSeq(0)))
        
    def blosum(self, col1, col2, d, e):
        blosum62 = Bio.SubsMat.MatrixInfo.blosum62
        blosum = 0
        
        for i in range(0, len(col1)):
            for j in range(0, len(col2)):
                if((col1[i], col2[j]) in blosum62):
                    blosum += blosum62[col1[i], col2[j]]
                elif((col2[j], col1[i]) in blosum62):
                    blosum += blosum62[col2[j], col1[i]]
                elif(col1[i] == '-'):
                    if(col2[j] == '+'):
                        blosum += -d
                    if(col2[j] != '-'):
                        blosum += -e
                elif(col2[j] == '-'):
                    if(col1[i] == '+'):
                        blosum += -d
                    if(col1[i] != '-'):
                        blosum += -e
                        
        return blosum / (len(col1)*len(col2))     
    
    def scoreIndexIsfrom(self, bloc, d, e):
        
        m = np.zeros((len(self.getSeqs()[0]) + 1, len(bloc.getSeqs()[0]) + 1)).astype(int)
        ix = np.zeros((len(self.getSeqs()[0]) + 1, len(bloc.getSeqs()[0]) + 1)).astype(int)
        iy = np.zeros((len(self.getSeqs()[0]) + 1, len(bloc.getSeqs()[0]) + 1)).astype(int)
        isfrom = np.zeros((len(self.getSeqs()[0]) + 1, len(bloc.getSeqs()[0]) + 1)).astype(int)
        
        m[1, 0] = -d
        m[0, 1] = -d
        ix[1, 0] = -d
        ix[0, 1] = -d
        iy[1, 0] = -d
        iy[0, 1] = -d
        isfrom[0, 1] = -1
        isfrom[1, 0] = 1
        
        for i in range(1, len(self.getSeqs()[0])):     #Calcul de la première ligne
            m[i+1, 0] = m[i, 0] - e
            ix[i+1, 0] = ix[i, 0] - e
            iy[i+1, 0] = iy[i, 0] - e
            isfrom[i+1, 0] = 1
            
        for j in range(1, len(bloc.getSeqs()[0])):     #Calcul de la première colonne
            m[0, j+1] = m[0, j] - e
            ix[0, j+1] = ix[0, j] - e
            iy[0, j+1] = iy[0, j] - e
            isfrom[0, j+1] = -1
        
        
        index, maxi = [0,0], m[0][0]
        
        for i in range(0, len(self.getSeqs()[0])):     #Relation de récurrence
            for j in range(0, len(bloc.getSeqs()[0])):
                
                coli = self.getCol(i)
                colj = bloc.getCol(j)
                
                blosum = self.blosum(coli, colj, d, e)
                
                #Calcul de m[i+1, j+1]
                score = max([iy[i, j], m[i, j], ix[i, j]]) + blosum
                m[i+1, j+1] = score
                
                if(score >= maxi):
                    maxi = score
                    index = [i+1, j+1]
                
                #Calcul de ix[i+1, j+1] et iy[i+1, j+1]
                ix[i+1, j+1] = max(m[i, j + 1] - d, ix[i, j + 1] - e)
                iy[i+1, j+1] = max(m[i + 1, j] - d, iy[i + 1, j] - e)
                
                #Calcul de isfrom[i+1, j+1]
                isfromij = np.argmax([iy[i+1, j+1], m[i+1, j+1], ix[i+1, j+1]]) - 1
                isfrom[i+1, j+1] = isfromij
                
        return maxi, index, isfrom
        
    def add(self, bloc, d, e):
        maxi, index, isfrom = self.scoreIndexIsfrom(bloc, d, e)
        
        seqs = []
        for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
            seqs += [Seq('')]
        
        while(0 not in index):
            if(isfrom.item(tuple(index)) == 0):
                index = [index[0] - 1, index[1] - 1]
                for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
                    if(i < self.getNbSeqs()):
                        seqs[i] = self.getSeq(i)[index[0]] + seqs[i]
                    else:
                        seqs[i] = bloc.getSeq(i - self.getNbSeqs())[index[1]] + seqs[i]
            
            elif(isfrom.item(tuple(index)) == 1):
                index = [index[0] - 1, index[1]]
                for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
                    if(i < self.getNbSeqs()):
                        seqs[i] = self.getSeq(i)[index[0]] + seqs[i]
                    else:
                        seqs[i] = '-' + seqs[i]
                
            elif(isfrom.item(tuple(index)) == -1):
                index = [index[0], index[1] - 1]
                for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
                    if(i < self.getNbSeqs()):
                        seqs[i] = '-' + seqs[i]
                    else:
                        seqs[i] = bloc.getSeq(i - self.getNbSeqs())[index[1]] + seqs[i]         
        
        self.seqs = seqs     
        self.nbSeqs = self.getNbSeqs() + bloc.getNbSeqs()
        self.score = maxi
        return maxi #retourner le score entre les deux
    
    def alignementscore(self, bloc, d, e):
        maxi, index, isfrom = self.scoreIndexIsfrom(bloc, d, e)
        return maxi #retourner le score entre les deux

#seq = Seq("WWAGCATTTGGCTGG")
#bloc1 = bloc(seq)
#bloc1.show()
#
#bloc2 = bloc(Seq("AGATGACTACCCT"))
#score = bloc1.alignementscore(bloc2, 10, 0.5)
#print(score)
#bloc1.show()
