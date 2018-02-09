# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:41:57 2018

@author: romai
"""
import Bio.SubsMat.MatrixInfo

class aminoAcidScorer:
    
    def __init__(self, methodName, params):
        self.name = methodName
        self.params = params
        
    def getName(self):
        return self.name
    
    def getParams(self):
        return self.params
    
    def blosum62(self, aa1, aa2):
        blosumMatrix = Bio.SubsMat.MatrixInfo.blosum62
        
        if((aa1["name"], aa2["name"]) in blosumMatrix):
            return blosumMatrix[aa1["name"], aa2["name"]]
        
        elif((aa2["name"], aa1["name"]) in blosumMatrix):
            return blosumMatrix[aa2["name"], aa1["name"]]
        
        elif(aa1["name"] == '-'):
            if(aa2["name"] == '+'):
                return -self.getParams()["openGap"]
            elif(aa2["name"] != '-'):
                return -self.getParams()["extendGap"]
            return 0
        
        elif(aa2["name"] == '-'):
            if(aa1["name"] == '+'):
                return -self.getParams()["openGap"]
            elif(aa1["name"] != '-'):
                return -self.getParams()["extendGap"]
            return 0
            
        elif(aa1["name"] == '+'):
            return -self.getParams()["openGap"]
        
        elif(aa2["name"] == '+'):
            return -self.getParams()["openGap"]
        
        print("WARNING - Can't compute the score for", aa1, aa2, "in 'blosum62' function")
        return 0

    #on met des multiplicateurs qui dilatent cette affaire
    def croisementPremier(self, aa1, aa2):
        result_temp = self.blosum62(aa1, aa2)

        #ma petite modif
        coef1 = 1 - aa1["enfouissement"] + 0.5
        coef2 = 1 - aa1["enfouissement"] + 0.5
        coefmixte = 3 * (aa1["struct"] != aa2["struct"]) + 1

        return coef1*coef2*coefmixte*result_temp
    
        print("WARNING - Can't compute the score for", aa1, aa2, "in 'blosum62' function")
        return 0
            
    def computeScore(self, aa1, aa2):
        if(self.getName() == "blosum62"):
            return self.blosum62(aa1, aa2)
        
        if (self.getName() == "blosum62mixte"):
            return self.croisementPremier(aa1, aa2)
        
        return 0