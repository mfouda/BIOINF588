# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:19:05 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo
from seqStruct import seqStruct

class bloc:
    
    def __init__(self, seq, name = 0):
        self.nbSeqs = 1
        self.seqs = [seq]
        self.dendo = ['']
        self.score = np.nan
        
    def getNbSeqs(self):
        return self.nbSeqs
    
    def getSeqs(self):
        return self.seqs
    
    def getSeq(self, i):
        return self.getSeqs()[i]
    
    def getNames(self):
        names = []
        for i in range(0, self.getNbSeqs()):
            names += [self.getSeq(i).getName()]
        return names
    
    def setNames(self, names):
        for i in range(0, self.getNbSeqs()):
            self.seqs[i].setName(names[i])
        
    def getName(self, i):
        return self.getSeq(i).getName()
    
    def getDendos(self):
        return self.dendo
    
    def getDendo(self, i):
        return self.getDendos()[i]
    
    def getCol(self, i):
        col = []
        for k in range(0, self.getNbSeqs()):
            col += [self.getSeq(k).getAminoAcid(i)]
        return col
    
    def show(self):
        print('#'*30)
        print('#######      BLOC      #######')
        print('#'*30)
        if(self.getNbSeqs() == 1):
            print('The bloc has 1 sequence')
        else:
            print('The bloc has '+str(self.getNbSeqs())+' sequences')
            
        print('The alignement is of lenght',len(self.getSeq(0)))
        if(not np.isnan(self.score)):
            print('The last merging score is',self.score)
        print('-'*30)
        if('aSequenceHasNoName' in self.getNames()):
            for i in range(0, self.getNbSeqs()):
                print(self.getSeq(i).printSeq(), self.getDendo(i))
        else:
             for i in range(0, self.getNbSeqs()):
                print(self.getSeq(i).printSeq(), self.getName(i), self.getDendo(i))   
        print('#'*30)
        
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
            seqs += [seqStruct()]
        
        minus, plus = dict(), dict()
        minus["name"], plus["name"] = "-", "+"
        minus["struct"], plus["struct"] = "", ""
        minus["enfouissement"], plus["enfouissement"] = 0, 0
                        
        while(0 not in index):
            if(isfrom.item(tuple(index)) == 0):
                index = [index[0] - 1, index[1] - 1]
                for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
                    if(i < self.getNbSeqs()):
                        seqs[i].addAminoAcidBefore(self.getSeq(i).getAminoAcid(index[0]))
                        #seqs[i] = self.getSeq(i)[index[0]] + seqs[i]
                    else:
                        seqs[i].addAminoAcidBefore(bloc.getSeq(i - self.getNbSeqs()).getAminoAcid(index[1]))
                        #seqs[i] = bloc.getSeq(i - self.getNbSeqs())[index[1]] + seqs[i]
            
            elif(isfrom.item(tuple(index)) == 1):
                index = [index[0] - 1, index[1]]
                for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
                    if(i < self.getNbSeqs()):
                        seqs[i].addAminoAcidBefore(self.getSeq(i).getAminoAcid(index[0]))
                        #seqs[i] = self.getSeq(i)[index[0]] + seqs[i]
                    else:
                        seqs[i].addAminoAcidBefore(minus)
                        #seqs[i] = '-' + seqs[i]
                
            elif(isfrom.item(tuple(index)) == -1):
                index = [index[0], index[1] - 1]
                for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
                    if(i < self.getNbSeqs()):
                        seqs[i].addAminoAcidBefore(minus)
                        #seqs[i] = '-' + seqs[i]
                    else:
                        seqs[i].addAminoAcidBefore(bloc.getSeq(i - self.getNbSeqs()).getAminoAcid(index[1]))
                        #seqs[i] = bloc.getSeq(i - self.getNbSeqs())[index[1]] + seqs[i]         
        
        for i in range(0, len(seqs)):
            if(seqs[i][0] == '-'):
                seqs[i].setAminoAcid(0, plus)
                #seqs[i] = '+' + seqs[i][1:]
            for k in range(0, len(seqs[0])-1):
                if(seqs[i].getAminoAcid(k+1)["name"] == '+' and (seqs[i].getAminoAcid(k)["name"] == '+' or seqs[i].getAminoAcid(k)["name"] == '-')):
                    seqs[i].setAminoAcid(k+1, minus)
                    #seqs[i] = seqs[i][:k+1] + '-' + seqs[i][k+2:]
                if(seqs[i].getAminoAcid(k+1)["name"] == '-' and (seqs[i].getAminoAcid(k)["name"] != '-' and seqs[i].getAminoAcid(k)["name"] != '+')):
                    seqs[i].setAminoAcid(k+1, minus)
                    #seqs[i] = seqs[i][:k+1] + '+' + seqs[i][k+2:]
        
        self.seqs = seqs        
        self.nbSeqs = self.getNbSeqs() + bloc.getNbSeqs()
        self.score = maxi
        
        
        dendo1 = self.getDendos()
        dendo2 = bloc.getDendos()
        ld1 = len(dendo1[0])
        ld2 = len(dendo2[0])        
        if(ld1 > ld2):
            for i in range(0, len(dendo1)):
                if(i == int(len(dendo1)/2)):
                    dendo1[i] += '-,'
                elif(i > int(len(dendo1)/2)):
                    dendo1[i] += ' |'
                else:
                    dendo1[i] += '  '
            for i in range(0, len(dendo2)):
                if(i == int(len(dendo2)/2)):
                    dendo2[i] = dendo2[i] + '-'*(ld1 - ld2) + "-'"
                elif(i > int(len(dendo2)/2)):
                    dendo2[i] = dendo2[i] + ' '*(ld1 - ld2) + '  '
                else:
                    dendo2[i] = dendo2[i] + ' '*(ld1 - ld2) + ' |'
        else:
            for i in range(0, len(dendo1)):
                if(i == int(len(dendo1)/2)):
                    dendo1[i] = dendo1[i] + '-'*(ld2 - ld1) + '-,'
                elif(i > int(len(dendo1)/2)):
                    dendo1[i] = dendo1[i] + ' '*(ld2 - ld1) + ' |'
                else:
                    dendo1[i] = dendo1[i] + ' '*(ld2 - ld1) + '  '
            for i in range(0, len(dendo2)):
                if(i == int(len(dendo2)/2)):
                    dendo2[i] += "-'"
                elif(i > int(len(dendo2)/2)):
                    dendo2[i] += '  '
                else:
                    dendo2[i] += ' |'
        self.dendo = dendo1 + dendo2
        self.setNames(self.getNames() + bloc.getNames())
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
#print(type(bloc1))
