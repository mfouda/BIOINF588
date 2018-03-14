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
import numpy as np
from seqStruct import seqStruct
from blocalignments import aligne_multiple
from score import aminoAcidScorer
import urllib.request as urllib
import os as os
from eval_SPS import SPS_romainTuned
import random
import time
import warnings
import pandas as pd
import datetime
warnings.filterwarnings("ignore")

from threading import Thread

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
                 
            for seq in SEQSmsffn:
                ID = []
                for k in range(0, seq.getLength()):
                    if(seq.getAminoAcid(k)["name"] != "-"):
                        ID += [seq.getAminoAcid(k)["id"]]
                if(not sum([ID[k] == k for k in range(0, len(ID))])):
                    print("WARNING - problème id des séquences", seq.getName(), ID)   
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
        #globalparams["openGap"] = [float(oGm.get()), float(oGp.get())]
        #globalparams["extendGap"] = [float(eGm.get()), float(eGp.get())]

        if "enf_mixte" not in globalparams:
            globalparams["enf_mixte"]= [0,0]

        if "vrac_mixte" not in globalparams:
            globalparams["vrac_mixte"] = [2,2]
        if "helice_mixte" not in globalparams:
            globalparams["helice_mixte"]= [4,4]
        if "struct_mixte" not in globalparams:
            globalparams["struct_mixte"] = [1,1]

        if "enf_propre" not in globalparams:
            globalparams["enf_propre"]= [0,0]

        if "struct_propre" not in globalparams:
            globalparams["struct_propre"] = [0,0]
        if "helice_propre" not in globalparams:
            globalparams["helice_propre"]= [4,4]

        if "openGap" not in globalparams:
            globalparams["openGap"] = [19,19]
        if "extendGap" not in globalparams:
            globalparams["extendGap"] = [2,2]

        print("Parameters : " + str(globalparams))
        
        PD = pd.DataFrame(columns = ["seqName", "SPS", "time", "iter"] + [k for k in globalparams.keys()], index=np.arange(0, int(niter.get()) * len(SEQS)))    

        SPS = 0
        bestparams = dict()
        
        ii = 0
        for i in range(0, int(niter.get())):
            params = dict()
            for k, v in globalparams.items():
                params[k] = round(10*(v[0] + (v[1] - v[0])*random.uniform(0, 1)))/10
            print("Iteration", i, str(params))
            
            if(multithreading.get()):
                keys = list(SEQS.keys())
                T = [ALIGNEMENT_SCORE([SEQSmsf[keys[k]].copy(), SEQS[keys[k]].copy(), aminoAcidScorer(str(sName.get()), params), keys[k]]) for k in range(0, len(keys))]
                [T[k].start() for k in range (0, len(keys))]
                [T[k].join() for k in range (0, len(keys))]
                tmpSPS = sum(T[k].getScore() for k in range(0, len(keys))) / len(SEQS)
            
            else:
                tmpSPS = 0
                scorer = aminoAcidScorer(str(sName.get()), params)
                for kk in SEQS.keys():
                    tt = time.time()
                    b = aligne_multiple(SEQS[kk], scorer)
                    tmptmpSPS = SPS_romainTuned(SEQSmsf[kk], b)
                    
                    PD.loc[ii]["seqName"] = kk
                    PD.loc[ii]["iter"] = i
                    PD.loc[ii]["SPS"] = tmptmpSPS
                    PD.loc[ii]["time"] = time.time() - tt
                    for k in params.keys():
                        PD.loc[ii][k] = params[k]
                    ii += 1
                    print(kk, tmptmpSPS, time.time() - tt)
                    tmpSPS += tmptmpSPS / len(SEQS)
            print("Iteration", i, "SPS =", int(100000*tmpSPS)/100000)
            if(tmpSPS > SPS):
                SPS = tmpSPS
                bestparams = params
          
        d = str(datetime.datetime.now())
        d = d[:4] + "_" + d[5:7] + "_" + d[8:10] + "_" + d[11:13] + "_" + d[14:16] + "_" + d[17:19]
        PD.to_csv("RESULT_" + d + ".csv", index = False)
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
    oGm.set(10)
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
    
    multithreading = tk.IntVar()
    c = tk.Checkbutton(FrameS, text="Multi threading", variable=multithreading, onvalue = True, offvalue = False, state=tk.DISABLED)
    c.pack()
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