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
        
        elif(aa1["name"] == '+'):
            return -self.getParams()["openGap"]
        
        elif(aa2["name"] == '+'):
            return -self.getParams()["openGap"]
        
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
            
        print("WARNING - Can't compute the score for", aa1, aa2, "in 'blosum62' function")
        return 0

    #on met des multiplicateurs qui dilatent cette affaire
    def blosum62mixte(self, aa1, aa2):
        result_temp = self.blosum62(aa1, aa2)

        #ma petite modif
        enf = self.getParams()["enf"]
        struct = self.getParams()["struct"]
        vrac = self.getParams()["vrac"]

        coef1_enf = 1.5 - aa1["enfouissement"]
        coef2_enf = 1.5 - aa2["enfouissement"]
        coef_enf = enf * coef1_enf * coef2_enf

        coef_mixte = 1
        if ((aa1["struct"] == aa2["struct"])) :
            if((aa1["struct"] != "V")):
                coef_mixte = vrac
            else :
                coef_mixte = struct

        return (1 + coef_enf + coef_mixte) * result_temp
    
    def testRomain(self, aa1, aa2):
        coefenf = 1 - (aa1["enfouissement"] - aa2["enfouissement"])**2
        enf = self.getParams()["enf"] * coefenf
        
        if(aa1["struct"] == aa2["struct"]):
            if(aa1["struct"] == "V"):
                coefmixte = self.getParams()["m1"]
            else:
                coefmixte = self.getParams()["m2"]
        else:
            coefmixte = 1
            
        coefmixte *= self.getParams()["mixte"]
        
        return enf + coefmixte + self.blosum62(aa1, aa2)
            
    def computeScore(self, aa1, aa2):
        if(self.getName() == "blosum62"):
            return self.blosum62(aa1, aa2)
        
        if (self.getName() == "blosum62mixte"):
            return self.blosum62mixte(aa1, aa2)
        
        if (self.getName() == "testRomain"):
            return self.testRomain(aa1, aa2)
        
        return 0