# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:41:57 2018

@author: romai
"""
import Bio.SubsMat.MatrixInfo

class aminoAcidScorer:
    
    def __init__(self, methodName, params):
        self.name = methodName

        if "enf_mixte" not in params:
            params["enf_mixte"]=1

        if "vrac_mixte" not in params:
            params["vrac_mixte"] = 2
        if "helice_mixte" not in params:
            params["helice_mixte"]= 4
        if "struct_mixte" not in params:
            params["struct_mixte"] = 1

        if "enf_propre" not in params:
            params["enf_propre"]= 1

        if "struct_propre" not in params:
            params["struct_propre"] = 1
        if "helice_propre" not in params:
            params["helice_propre"]= 4

        if "openGap" not in params:
                params["openGap"] = 19
        if "extendGap" not in params:
                params["extendGap"] = 2
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

    #Ici on définit des fonctions possibles pour les fonctions d'évaluation, propres et mixtes, basées sur l'enfouissement et la structure
    def enf_propre(self, aa1, aa2, enf):
        coef1_enf = 1.5 - aa1["enfouissement"]
        coef2_enf = 1.5 - aa2["enfouissement"]
        coef_enf = enf * (coef1_enf + coef2_enf)
        return coef_enf

    def enf_mixte(self, aa1, aa2, enf):
        coef_enf = 1 - (aa1["enfouissement"] - aa2["enfouissement"]) ** 2
        coef_enf *= enf
        return coef_enf

    def struct_propre(self, aa1, aa2, struct, helice):
        coef1 = 1
        if (aa1["struct"] != "V"):
            coef1 = helice
            
        coef2 = 1
        if (aa2["struct"] != "V"):
            coef2 = helice
        return struct * (coef1 + coef2)

    def struct_mixte(self, aa1, aa2, vrac, helice, struct):
        if (aa1["struct"] == aa2["struct"]):
            if (aa1["struct"] == "V"):
                coef_struct = vrac
            else:
                coef_struct = helice
        else:
            coef_struct= 1

        coef_struct *= struct
        return coef_struct

    #On définit des scorers
    #scorer purement mixte
    def coef_mixte(self, aa1, aa2):
        
        if(aa1["name"] in ["-", "+"] or aa2["name"] in ["-", "+"]):
            return 0
        else:
            enf = self.getParams()["enf_mixte"]
    
            vrac = self.getParams()["vrac_mixte"]
            helice = self.getParams()["helice_mixte"]
            struct = self.getParams()["struct_mixte"]
    
            #terme mixte d'enfouissement
            coef_enf = self.enf_mixte(aa1, aa2, enf)
    
            #terme mixte structural
            coef_struct = self.struct_mixte(aa1, aa2, vrac, helice, struct)
    
            return coef_enf + coef_struct

    def mixte(self, aa1, aa2):
        return self.blosum62(aa1, aa2) + self.coef_mixte(aa1, aa2)

    #scorer purement propre
    def coef_propre(self, aa1, aa2):
        
        if(aa1["name"] in ["-", "+"] or aa2["name"] in ["-", "+"]):
            return 0
        else:
            enf = self.getParams()["enf_propre"]
    
            struct = self.getParams()["struct_propre"]
            helice = self.getParams()["helice_propre"]
    
            #terme propre d'enfouissement
            coef_enf = self.enf_propre(aa1, aa2, enf)
    
            # terme propre structural
            coef_struct = self.struct_propre(aa1, aa2, struct, helice)
    
            return coef_enf + coef_struct

    def propre(self, aa1, aa2):
        return self.blosum62(aa1, aa2) + self.coef_propre(aa1, aa2)

    #scorer avec effets mixte et propre
    def total(self, aa1, aa2):
        return self.blosum62(aa1, aa2) + self.coef_propre(aa1, aa2) + self.coef_mixte(aa1, aa2)


    def testRomain(self, aa1, aa2):
        coefenf = 1 - (aa1["enfouissement"] - aa2["enfouissement"])**2
        enf = self.getParams()["enf_mixte"] * coefenf
        
        if(aa1["struct"] == aa2["struct"]):
            if(aa1["struct"] == "V"):
                coefmixte = self.getParams()["vrac_mixte"]
            else:
                coefmixte = self.getParams()["helice_mixte"]
        else:
            coefmixte = 1
            
        coefmixte *= self.getParams()["struct_mixte"]
        
        return enf + coefmixte + self.blosum62(aa1, aa2)
    
    def score_Jas(self, aa1, aa2) :
        enf1 = aa1["enfouissement"]
        enf2 = aa2["enfouissement"]
        enf = self.getParams()["enf"]*(0.5 + enf1)*(0.5 + enf2)*(enf1 - enf2)**2
        
        struct = self.getParams()["struct"]
        s1 = aa1["struct"]
        s2 = aa2["struct"]
        if s1 == 'V' :
            if s2 == 'V' :
                struct *= self.getParams()["p1"]
            else :
                struct *= self.getParams()["p4"]
        else :
            if s1 == s2 :
                struct *= self.getParams()["p2"]
            else :
                struct *= self.getParams()["p3"]
        #a priori p1<p2<p3<p4
        
        return (1 + enf + struct)*self.blosum62(aa1, aa2)
    
    
    def computeScore(self, aa1, aa2):
        if(self.getName() == "blosum62"):
            return self.blosum62(aa1, aa2)
        if (self.getName() == "mixte"):
            return self.mixte(aa1, aa2)
        if (self.getName() == "propre"):
            return self.propre(aa1, aa2)
        if (self.getName() == "total"):
            return self.total(aa1, aa2)
        if (self.getName() == "testRomain"):
            return self.testRomain(aa1, aa2)
        if (self.getName() == "score_Jas"):
            return self.score_Jas(aa1, aa2)
        return 0