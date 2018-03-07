# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 10:46:13 2018

@author: romai
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:19:10 2018

@author: romai
"""

# On importe Tkinter
import tkinter as tk
import tkinter.filedialog as tkfd
#import random as random
from seqStruct import seqStruct
from blocalignments import aligne_multiple
from score import aminoAcidScorer
import urllib.request as urllib
import os as os
from eval_SPS import SPS_romainTuned
import random

import warnings
warnings.filterwarnings("ignore")

import sys
from threading import Thread
import time

class Afficheur(Thread):

    """Thread chargé simplement d'afficher une lettre dans la console."""

    def __init__(self, lettre):
        Thread.__init__(self)
        self.lettre = lettre
        self.value = 0

    def run(self):
        """Code à exécuter pendant l'exécution du thread."""
        self.value = self.lettre ** 2
        attente = 0.2
        attente += random.randint(1, 120) / 100
        time.sleep(attente)
        print(self.value)

    def getValue(self):
        return self.value

class ALIGNEMENT_SCORE(Thread):

    """Thread chargé simplement d'afficher une lettre dans la console."""

    def __init__(self, VAL):
        Thread.__init__(self)
        self.msf = VAL[0]
        self.seqs = VAL[1]
        self.scorer = VAL[2]
        self.label = VAL[3]
        self.score = 0

    def run(self):
        """Code à exécuter pendant l'exécution du thread."""
        self.score = SPS_romainTuned(self.msf, aligne_multiple(self.seqs, self.scorer))
        print(self.label, int(1000*self.score)/1000)
        
    def getScore(self):
        return self.score

def showProtein(name):
    seq = seqStruct("pdb/" + name + ".pdb")
    
    # Création de la fenêtre principale (main window)
    Mafenetre = tk.Tk()
    Mafenetre.title('Protein '+name)   
    
    # Création d'un widget Canvas (zone graphique)
    Largeur = 10*len(seq.toString())
    Hauteur = 100
    Canevas = tk.Canvas(Mafenetre, width = Largeur, height =Hauteur, bg ='white')
    Canevas.pack(padx =5, pady =5)
    
    stringseq = seq.toString()
    for i in range(0, len(stringseq)):
        aa = seq.getAminoAcid(i)
        if(aa["struct"] == "H"):
            fill = "red"
        elif(aa["struct"] == "F"):
            fill = "blue"
        else:
            fill = "black"
        Canevas.create_text(10 + (i / (len(stringseq) - 1)) * (Largeur - 20), 15, 
                            text = stringseq[i], font=("Helvetica", 10), fill = fill)
        if(i != len(stringseq) - 1):
            Canevas.create_line(10 + (i / (len(stringseq) - 1)) * (Largeur - 20), 95 - 75 * seq.getAminoAcid(i)["enfouissement"],
                                10 + ((i + 1) / (len(stringseq) - 1)) * (Largeur - 20), 95 - 75 * seq.getAminoAcid(i + 1)["enfouissement"])
    # Création d'un widget Button (bouton Quitter)
    BoutonQuitter = tk.Button(Mafenetre, text ='Quitter', command = Mafenetre.destroy)
    BoutonQuitter.pack(side = tk.LEFT, padx = 5, pady = 5)
    
    Mafenetre.mainloop()  
    
def launchInterface():
    
    def Reset():
        global SEQS, params
        """ Efface la zone graphique et réinitialise tout"""
        print("RESETING all variables")
        SEQS = []
        params = dict()
        ButtonA['state'] = 'disable'
        ButtonAps['state'] = 'disable'
    
    def msfToSeqs():
        global filename, SEQSmsf
        SEQSmsf = dict()
        for fn in filename:
            fn = fn[:-4] + ".msf"
            print("Ouverture de " + fn[-11:])
            lines = open(fn, 'r').readlines()
            
            names = []
            i = 6
            while(lines[i][:5] == ' Name'):
                names += [lines[i][7:11]]
                i += 1
            
            SEQSmsffn = []
            last = dict()
            for i in range(0, len(names)):
                SEQSmsffn += [seqStruct()]
                SEQSmsffn[-1].setName(names[i])
                last[names[i]] = 0
            
            for line in lines:
                if(line[:4] in names):
                    index = names.index(line[:4])
                    words = line.split()
                    for word in words[1:] :
                        for aa in word :
                            d = dict()
                            d["name"] = aa
                            d["struct"] = "V"
                            if(aa != "."):
                                d["id"]=last[line[:4]]
                                last[line[:4]] += 1
                            else:
                                d["name"] = "-"
                            SEQSmsffn[index].addAminoAcidAfter(d)
                            
#            for i in range(0, len(SEQSmsffn)):
#                seq = seqStruct("pdb/" + SEQSmsffn[i].getName() + ".pdb")
#                n = 0
#                for k in range(0, SEQSmsffn[i].getLength()):
#                    if("id" in SEQSmsffn[i].getAminoAcid(k)):
#                        aa = SEQSmsffn[i].getAminoAcid(k).copy()
#                        aaseq = seq.getAminoAcid(n).copy()
#                        aa["struct"] = aaseq["struct"]
#                        aa["enfouissement"] = aaseq["enfouissement"]
#                        SEQSmsffn[i].setAminoAcid(k, aa)
#                        n += 1
            SEQSmsf[fn[-11:-4]] = SEQSmsffn            

    def searchfilename():
        global filename
        filename = [tkfd.askopenfilename(initialdir = "../RV11/", title="Ouvrir un .tfa", filetypes=[('tfa files','.tfa'),('all files','.*')])]
        useTFA()
        
    def searchdirectory():
        global filename
        filedir = tkfd.askdirectory(initialdir = "../", title="Ouvrir un dossier")
        filename = []
        for fn in os.listdir(filedir):
            if(fn[-4:] == ".tfa"):
                filename += [filedir + "/" + fn]
        if(len(filename) != 0):
            useTFA()
        else:
            filename = []
            print("Pas de .tfa dans le dossier")
        
    def useTFA():
        global SEQS, globalparams, filename
        
        SEQS = dict()
        for i in range(0, len(filename)):
            fn = filename[i]
            try:
                SEQSfn = []
                print("Recherche de séquences dans " + fn[-11:])
                lines = open(fn, 'r').readlines()
                #lines = [lines[i][:-2] for i in range(0, len(lines))]
                for line in lines:
                    if(line[0] == ">"):
                        if(line[1:5] + ".pdb" not in os.listdir("pdb/")):
                            print(" - - - - > Downloading " + line[1:5] + ".pdb ...")
                            path = "https://files.rcsb.org/download/" + line[1:5] + ".pdb"
                            urllib.urlretrieve(path, "pdb/" + line[1:5] + ".pdb")
                        SEQSfn += [seqStruct("pdb/" + line[1:5] + ".pdb")]
                SEQS[fn[-11:-4]] = SEQSfn
            except (KeyError):
                print(fn[-11:] + " ERROR - Fichier .tfa illisible ")
                filename = filename[:i] + filename[i+1:]
            except (urllib.HTTPError):
                print(fn[-11:] + " ERROR - Fichier .pdb impossible a télécharger")
                filename = filename[:i] + filename[i+1:]
        globalparams = dict()
        ButtonA['state'] = 'normal'
        ButtonAps['state'] = 'normal'

    def addParams():
        global globalparams
        globalparams[str(key.get())] = [float(valuem.get()), float(valuep.get())]
        print("Adding entry in params  {" + str(key.get()) + " : [" + str(valuem.get())+", "+str(valuep.get())+"]}")
        
    def optimize():
        global SEQS, globalparams, bloc, SEQSmsf, filename
        msfToSeqs()
        globalparams["openGap"] = [float(oGm.get()), float(oGp.get())]
        globalparams["extendGap"] = [float(oGm.get()), float(eGp.get())]
        
        print("Parameters : " + str(globalparams))
        
        SPS = 0
        bestparams = dict()
        
        for i in range(0, int(niter.get())):
            params = dict()
            for k, v in globalparams.items():
                params[k] = round(10*(v[0] + (v[1] - v[0])*random.uniform(0, 1)))/10
            print("Iteration", i, str(params))
            tmpSPS = 0
            T = []
            keys = [1, 2, 3, 4, 5, 6]
            
#            for k in range(0, len(keys)):
#                T += [Afficheur(keys[k])]
#            for k in range(0, len(keys)):
#                T[k].start()
#            for k in range(0, len(keys)):
#                T[k].join()
            
            scorer = aminoAcidScorer(str(sName.get()), params)
            keys = list(SEQS.keys())
            for k in range(0, len(keys)):
                key = keys[k]
                T += [ALIGNEMENT_SCORE([SEQSmsf[key], SEQS[key], scorer, key])]
            for k in range(0, len(keys)):
                T[k].start()
            for k in range(0, len(keys)):
                T[k].join()   
            for k in range(0, len(keys)):
                tmpSPS += T[k].getScore() / len(SEQS)
#            for k in SEQS.keys():
#                tmptmpSPS = SPS_romainTuned(SEQSmsf[k], aligne_multiple(SEQS[k], scorer))
#                tmpSPS += tmptmpSPS / len(SEQS)
#                print(i, k, int(1000*tmptmpSPS)/1000)
            print("Iteration", i, "SPS =", int(100000*tmpSPS)/100000)
            if(tmpSPS > SPS):
                SPS = tmpSPS
                bestparams = params
                
        print(bestparams, SPS)
            
    # Création de la fenêtre principale (main window)
    Mafenetre = tk.Tk()
    Mafenetre.title('Protein alignment')
    
    #Widget Protein
    FramePName = tk.Frame(Mafenetre,borderwidth=2,relief=tk.GROOVE)
    FramePName.pack(side = tk.LEFT, padx=2,pady=2)
    LabelP = tk.Label(FramePName,text="Fichiers à analyser")
    LabelP.pack(padx=2,pady=2)
    
    ButtonTFA = tk.Button(FramePName,text="Utiliser 1 .tfa",fg='navy',command=searchfilename)
    ButtonTFA.pack(padx=2,pady=2)
    ButtonTFAm = tk.Button(FramePName,text="Utiliser plusieurs .tfa",fg='navy',command=searchdirectory)
    ButtonTFAm.pack(padx=2,pady=2)
    
    #Widget Scorer
    FrameS = tk.Frame(Mafenetre,borderwidth=2,relief=tk.GROOVE)
    FrameS.pack(side = tk.LEFT, padx=2,pady=2)
    LabelS= tk.Label(FrameS,text="Scorer")
    LabelS.pack(side = tk.TOP, padx=2,pady=2)
    sName = tk.StringVar()
    sName.set("blosum62")
    TextVarsName = tk.Entry(FrameS, textvariable=sName, width=20)
    TextVarsName.pack(side = tk.TOP, padx=2,pady=2)
    
    FrameSS = tk.Frame(FrameS,borderwidth=1)
    FrameSS.pack(padx=2,pady=2)
    oGp = tk.StringVar()
    oGp.set(25)
    TextVaroGp = tk.Entry(FrameSS, textvariable=oGp, width=10)
    oGm = tk.StringVar()
    oGm.set(2)
    TextVaroGm = tk.Entry(FrameSS, textvariable=oGm, width=10)
    #TextVaroG.pack(padx=2,pady=2)
    TextVaroGp.grid(row =2, column =0, padx=1,pady=1)
    TextVaroGm.grid(row =1, column =0, padx=1,pady=1)

    eGp = tk.StringVar()
    eGp.set(5)
    TextVareGp = tk.Entry(FrameSS, textvariable=eGp, width=10)
    eGm = tk.StringVar()
    eGm.set(0.5)
    TextVareGm = tk.Entry(FrameSS, textvariable=eGm, width=10)
    #TextVaroG.pack(padx=2,pady=2)
    TextVareGp.grid(row =2, column =1, padx=1,pady=1)
    TextVareGm.grid(row =1, column =1, padx=1,pady=1)
    #TextVareG.pack(padx=2,pady=2)
    #TextVareG.grid(row =1, column =1, padx=1,pady=1)
    tk.Label(FrameSS, text = 'Open gap').grid(row =0, column =0, padx=1,pady=1)
    tk.Label(FrameSS, text = 'Extend gap').grid(row =0, column =1, padx=1,pady=1)
    
    FrameAps = tk.Frame(FrameS,borderwidth=2,relief=tk.GROOVE)
    FrameAps.pack(padx=2,pady=2)
    
    FrameApss = tk.Frame(FrameAps,borderwidth=2)
    FrameApss.pack(padx=2,pady=2)
    key = tk.StringVar()
    key.set("key")
    TextVarkey = tk.Entry(FrameApss, textvariable=key, width=10)
    TextVarkey.grid(row =0, column =0, padx=1,pady=1)
    valuep = tk.StringVar()
    valuep.set(1)
    TextVarvaluep = tk.Entry(FrameApss, textvariable=valuep, width=10)
    TextVarvaluep.grid(row =0, column =1, padx=1,pady=1)
    valuem = tk.StringVar()
    valuem.set(2)
    TextVarvaluem = tk.Entry(FrameApss, textvariable=valuem, width=10)
    TextVarvaluem.grid(row =0, column =2, padx=1,pady=1)
    ButtonAps = tk.Button(FrameAps,text="Add parameter",fg='navy',command=addParams, state=tk.DISABLED)
    ButtonAps.pack(padx=2,pady=2)
    
    FrameSO = tk.Frame(FrameS,borderwidth=2)
    FrameSO.pack(padx=2,pady=2)
    tk.Label(FrameSO, text = 'iterations').grid(row =0, column =0, padx=1,pady=1)
    niter = tk.StringVar()
    niter.set(5)
    TextVarvalueniter = tk.Entry(FrameSO, textvariable=niter, width=10)
    TextVarvalueniter.grid(row =0, column =1, padx=1,pady=1)
    ButtonA = tk.Button(FrameS,text="Optimiser",fg='navy',command=optimize, state=tk.DISABLED)
    ButtonA.pack(side = tk.BOTTOM, padx=2,pady=2)
    
    FrameABU = tk.Frame(Mafenetre,borderwidth=2)
    FrameABU.pack(side = tk.LEFT, padx=2,pady=2)

    # Création d'un widget Button (bouton Effacer)
    BoutonReset = tk.Button(FrameABU, text ='Reset', command = Reset)
    BoutonReset.pack(padx = 5, pady = 5)
    # Création d'un widget Button (bouton Quitter)
    BoutonQuitter = tk.Button(FrameABU, text ='Quitter', command = Mafenetre.destroy)
    BoutonQuitter.pack(padx = 5, pady = 5)
    
    Mafenetre.mainloop()
    
launchInterface()