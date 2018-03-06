# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:19:05 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo
from seqStruct import seqStruct
from score import aminoAcidScorer
import pickle as pk
import random
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
 

class bloc:
    
    def __init__(self, seq):
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
    
    def getLength(self):
        return self.getSeq(0).getLength()
    
    def show(self):
        print('#'*30)
        print('#######      BLOC      #######')
        print('#'*30)
        if(self.getNbSeqs() == 1):
            print('The bloc has 1 sequence')
        else:
            print('The bloc has '+str(self.getNbSeqs())+' sequences')
            
        print('The alignement is of lenght', self.getLength())
        if(not np.isnan(self.score)):
            print('The last merging score is',self.score)
        print('-'*30)
        if('aSequenceHasNoName' in self.getNames()):
            for i in range(0, self.getNbSeqs()):
                print(self.getSeq(i).toString(), self.getDendo(i))
        else:
             for i in range(0, self.getNbSeqs()):
                print(self.getSeq(i).toString(), self.getName(i) + " "*(max([len(n) for n in self.getNames()]) - len(self.getName(i))), self.getDendo(i))   
        print('#'*30)
        
    def aminoAcidScore(self, col1, col2, scorer):
        aaScore = 0
        
        for i in range(0, len(col1)):
            for j in range(0, len(col2)):
                aaScore += scorer.computeScore(col1[i], col2[j])
                
        return aaScore / (len(col1)*len(col2))     
    
    def scoreIndexIsfrom(self, bloc, scorer):
        
        size = (self.getSeq(0).getLength() + 1, bloc.getSeq(0).getLength() + 1)
        m = np.zeros(size).astype(int)
        ix = np.zeros(size).astype(int)
        iy = np.zeros(size).astype(int)
        isfrom = np.zeros(size).astype(int)
        del size
        
        minus = dict({"name" : "-", "struct" : "V", "enfouissement" : 0})
        plus = dict({"name" : "+", "struct" : "V", "enfouissement" : 0})
        
        ds = self.aminoAcidScore(self.getCol(0), [plus], scorer)
        db = self.aminoAcidScore(bloc.getCol(0), [plus], scorer)
        m[1, 0] = ds
        m[0, 1] = db
        ix[1, 0] = ds
        ix[0, 1] = db
        iy[1, 0] = ds
        iy[0, 1] = db
        isfrom[0, 1] = -1
        isfrom[1, 0] = 1
        
        for i in range(1, self.getSeq(0).getLength()):     #Calcul de la première ligne
            e = self.aminoAcidScore(self.getCol(i), [minus], scorer)
            m[i+1, 0] = m[i, 0] + e
            ix[i+1, 0] = ix[i, 0] + e
            iy[i+1, 0] = iy[i, 0] + e
            isfrom[i+1, 0] = 1
            
        for j in range(1, bloc.getSeq(0).getLength()):     #Calcul de la première colonne
            e = self.aminoAcidScore(bloc.getCol(j), [minus], scorer)
            m[0, j+1] = m[0, j] + e
            ix[0, j+1] = ix[0, j] + e
            iy[0, j+1] = iy[0, j] + e
            isfrom[0, j+1] = -1
        
        #Inutile pour le moment
        index, maxi = [0,0], m[0][0]
        
        for i in range(0, self.getSeq(0).getLength()):     #Relation de récurrence
            coli = self.getCol(i)
            for j in range(0, bloc.getSeq(0).getLength()):
                colj = bloc.getCol(j)
                
                blosum = self.aminoAcidScore(coli, colj, scorer)
                
                #Calcul de m[i+1, j+1]
                score = max([iy[i, j], m[i, j], ix[i, j]]) + blosum
                m[i+1, j+1] = score
                
                if(score >= maxi and 
                   ((j + 1 == bloc.getSeq(0).getLength() and i != 0) or 
                    (i + 1 == self.getSeq(0).getLength() and j != 0))):
                    maxi = score
                    index = [i+1, j+1]
                
                #Calcul de ix[i+1, j+1] et iy[i+1, j+1]
                #ix[i+1, j+1] = max(m[i, j + 1] - scorer.getParams()["openGap"], ix[i, j + 1] - scorer.getParams()["extendGap"])
                #iy[i+1, j+1] = max(m[i + 1, j] - scorer.getParams()["openGap"], iy[i + 1, j] - scorer.getParams()["extendGap"])
                ix[i+1, j+1] = max(m[i, j + 1] + self.aminoAcidScore(coli, [plus], scorer),
                                  ix[i, j + 1] + self.aminoAcidScore(coli, [minus], scorer))
                iy[i+1, j+1] = max(m[i + 1, j] + self.aminoAcidScore(colj, [plus], scorer),
                                  iy[i + 1, j] + self.aminoAcidScore(colj, [minus], scorer))
                
                #Calcul de isfrom[i+1, j+1]
                isfrom[i+1, j+1] = np.argmax([iy[i+1, j+1], m[i+1, j+1], ix[i+1, j+1]]) - 1
            
        writer = ExcelWriter('pickleObjects/test' + str(random.randint(0, 10000)) + '.xlsx')
        pd.DataFrame(m).to_excel(writer,'Sheet1',index=False)
        pd.DataFrame(ix).to_excel(writer,'Sheet2',index=False)
        pd.DataFrame(iy).to_excel(writer,'Sheet3',index=False)
        pd.DataFrame(isfrom).to_excel(writer,'Sheet4',index=False)
        writer.save()

        return maxi, index, isfrom
        
    def add(self, bloc, scorer):
        maxi, index, isfrom = self.scoreIndexIsfrom(bloc, scorer)
#        index = list(isfrom.shape)
#        index[0] -= 1
#        index[1] -= 1
        minus = dict({"name" : "-", "struct" : "", "enfouissement" : 0})
        plus = dict({"name" : "+", "struct" : "", "enfouissement" : 0})
        
        seqs = []
        for i in range(0, self.getNbSeqs() + bloc.getNbSeqs()):
            seqs += [seqStruct()]
            if(i < self.getNbSeqs()):
                seqs[-1].setName(self.getSeq(i).getName())
                for k in range(index[0] + 1, self.getSeq(0).getLength()):
                    seqs[-1].addAminoAcidAfter(self.getSeq(i).getAminoAcid(k))
                for k in range(index[1] + 1, bloc.getSeq(0).getLength()):
                    seqs[-1].addAminoAcidAfter(minus)
                
            else:
                seqs[-1].setName(bloc.getSeq(i - self.getNbSeqs()).getName())
                for k in range(index[0] + 1, self.getSeq(0).getLength()):
                    seqs[-1].addAminoAcidAfter(minus)
                for k in range(index[1] + 1, bloc.getSeq(0).getLength()):
                    seqs[-1].addAminoAcidAfter(bloc.getSeq(i - self.getNbSeqs()).getAminoAcid(k))
        
        minus = dict({"name" : "-", "struct" : "", "enfouissement" : 0})
        plus = dict({"name" : "+", "struct" : "", "enfouissement" : 0})
        
        #Rajouter la suite des séquences
        

        #while(0 not in index):
        while(index != [0, 0]):
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
#            if(seqs[i].getAminoAcid(0)["name"] == '-'):
#                seqs[i].setAminoAcid(0, plus)
                #seqs[i] = '+' + seqs[i][1:]
            for k in range(0, seqs[i].getLength()-1):
                if(seqs[i].getAminoAcid(k+1)["name"] == '+' and (seqs[i].getAminoAcid(k)["name"] == '+' or seqs[i].getAminoAcid(k)["name"] == '-')):
                    seqs[i].setAminoAcid(k+1, minus)
                    #seqs[i] = seqs[i][:k+1] + '-' + seqs[i][k+2:]
                if(seqs[i].getAminoAcid(k+1)["name"] == '-' and (seqs[i].getAminoAcid(k)["name"] != '-' and seqs[i].getAminoAcid(k)["name"] != '+')):
                    seqs[i].setAminoAcid(k+1, plus)
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
    
    def alignementscore(self, bloc, scorer):
        maxi, index, isfrom = self.scoreIndexIsfrom(bloc, scorer)
        return maxi #retourner le score entre les deux