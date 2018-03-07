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
import datetime
from seqStruct import seqStruct
from blocalignments import aligne_multiple
from score import aminoAcidScorer
import urllib.request as urllib
import os as os
from eval_SPS import SPS_romainTuned
import random

import warnings
warnings.filterwarnings("ignore")
    
def showBloc(bloc):
    # Création de la fenêtre principale (main window)
    Mafenetre = tk.Tk()
    
    LINES = []
    LINES += ['#The bloc has '+str(bloc.getNbSeqs())+' sequence'+"s"*(bloc.getNbSeqs() != 1)]
    LINES += ['#The alignement is of lenght ' + str(bloc.getSeq(0).getLength())]
    if(not np.isnan(bloc.score)):
        LINES += ['#The last merging score is ' + str(bloc.score)]
    LINES += ["#" + ''*30]
    if('aSequenceHasNoName' in bloc.getNames()):
        for i in range(0, bloc.getNbSeqs()):
            LINES += [bloc.getSeq(i).toString() + " " + bloc.getDendo(i)]
    else:
         for i in range(0, bloc.getNbSeqs()):
            LINES += [bloc.getSeq(i).toString() + " " + bloc.getName(i) + " "*(max([len(n) for n in bloc.getNames()]) - len(bloc.getName(i))) + " " + bloc.getDendo(i)]  
    
    # Création d'un widget Canvas (zone graphique)
    lenbloc = len(bloc.getSeq(0).toString())
    maxlen = max(len(str(LINES[i])) for i in range(0, len(LINES)))
    Largeur = 10*maxlen
    Hauteur = 20*len(LINES)
    Canevas = tk.Canvas(Mafenetre, width = Largeur, height =Hauteur, bg ='white')
    Canevas.pack(padx =5, pady =5)
    numbloc = 0
    for k in range(0, len(LINES)):
        if(LINES[k][0] == "#"):
           Canevas.create_text(10, 10 + (k / (len(LINES) - 1)) * (Hauteur - 20), anchor = tk.NW,
                                text = LINES[k][1:], font=("Helvetica", 10), fill = "black")
        else:
            for i in range(0, len(LINES[k])):
                if(i < lenbloc):
                    aa = bloc.getSeq(numbloc).getAminoAcid(i)
                    if(aa["struct"] == "H"):
                        fill = "red"
                    elif(aa["struct"] == "F"):
                        fill = "blue"
                    else:
                        fill = "black"
                else:
                    fill = "black"
                Canevas.create_text(10 + (i / (maxlen - 1)) * (Largeur - 20), 10 + (k / (len(LINES) - 1)) * (Hauteur - 20),
                                    text = LINES[k][i], font=("Helvetica", 10), fill = fill)
            numbloc += 1
            
    # Création d'un widget Button (bouton Quitter)
    BoutonQuitter = tk.Button(Mafenetre, text ='Quitter', command = Mafenetre.destroy)
    BoutonQuitter.pack(side = tk.LEFT, padx = 5, pady = 5)
    
    Mafenetre.mainloop()

def showSEQS(SEQS):
    # Création de la fenêtre principale (main window)
    Mafenetre = tk.Tk()
    
    # Création d'un widget Canvas (zone graphique)
    Largeur = 10*len(SEQS[0].toString())
    Hauteur = 20*len(SEQS)
    Canevas = tk.Canvas(Mafenetre, width = Largeur, height =Hauteur, bg ='white')
    Canevas.pack(padx =5, pady =5)
    
    maxlen = max(len(SEQS[k].toString()) for k in range(0, len(SEQS)))
    for k in range(0, len(SEQS)):
        stringseq = SEQS[k].toString()
        for i in range(0, len(stringseq)):
            aa = SEQS[k].getAminoAcid(i)
            if(aa["struct"] == "H"):
                fill = "red"
            elif(aa["struct"] == "F"):
                fill = "blue"
            else:
                fill = "black"
            Canevas.create_text(10 + (i / (maxlen - 1)) * (Largeur - 20), 10 + (k / (len(SEQS) - 1)) * (Hauteur - 20), 
                                text = stringseq[i], font=("Helvetica", 10), fill = fill)
            
    # Création d'un widget Button (bouton Quitter)
    BoutonQuitter = tk.Button(Mafenetre, text ='Quitter', command = Mafenetre.destroy)
    BoutonQuitter.pack(side = tk.LEFT, padx = 5, pady = 5)
    
    Mafenetre.mainloop()  
    
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
    
    logs = []
    
    def Reset():
        global SEQS, params
        """ Efface la zone graphique et réinitialise tout"""
        log(logs, "RESETING all variables")
        SEQS = []
        params = dict()
        ButtonA['state'] = 'disable'
        ButtonAps['state'] = 'disable'
        
    def Effacer():
        """ Efface la zone graphique """
        logger.delete(tk.ALL)
    
    def msfToSeqs():
        global filename, SEQSmsf
        SEQSmsf = dict()
        for fn in filename:
            fn = fn[:-4] + ".msf"
            log(logs, "Ouverture d'un .msf ...")
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
                            
            for i in range(0, len(SEQSmsffn)):
                seq = seqStruct("pdb/" + SEQSmsffn[i].getName() + ".pdb")
                n = 0
                for k in range(0, SEQSmsffn[i].getLength()):
                    if("id" in SEQSmsffn[i].getAminoAcid(k)):
                        aa = SEQSmsffn[i].getAminoAcid(k).copy()
                        aaseq = seq.getAminoAcid(n).copy()
                        aa["struct"] = aaseq["struct"]
                        aa["enfouissement"] = aaseq["enfouissement"]
                        SEQSmsffn[i].setAminoAcid(k, aa)
                        n += 1
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
            log(logs, "Pas de .tfa dans le dossier")
        
    def useTFA():
        global SEQS, globalparams, filename
        
        SEQS = dict()
        for i in range(0, len(filename)):
            fn = filename[i]
            try:
                SEQSfn = []
                log(logs, "Recherche de séquences dans le .tfa")
                lines = open(fn, 'r').readlines()
                #lines = [lines[i][:-2] for i in range(0, len(lines))]
                for line in lines:
                    if(line[0] == ">"):
                        if(line[1:5] + ".pdb" not in os.listdir("pdb/")):
                            log(logs, "Downloading " + line[1:5] + ".pdb ...")
                            path = "https://files.rcsb.org/download/" + line[1:5] + ".pdb"
                            urllib.urlretrieve(path, "pdb/" + line[1:5] + ".pdb")
                        else:
                            log(logs, "Protein " + line[1:5] + ".pdb already downloaded before")
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
            
    def log(logs, s):
        Effacer()
        logs += ["[" + str(datetime.datetime.now())[:-7] + "]    " + s]
        nrows = int(Hauteur / 10)
        for i in range(0, min(nrows, len(logs))):
            logger.create_text(5, 5 + 1.1*(i / nrows) * (Hauteur - 10), text = logs[-i-1], anchor = tk.NW, font=("Helvetica", 9))
    
    def addParams():
        global globalparams
        globalparams[str(key.get())] = [float(valuem.get()), float(valuep.get())]
        log(logs, "Adding entry in params  {" + str(key.get()) + " : [" + str(valuem.get())+", "+str(valuep.get())+"]}")
        
    def optimize():
        global SEQS, globalparams, bloc, SEQSmsf, filename
        msfToSeqs()
        globalparams["openGap"] = [float(oGm.get()), float(oGp.get())]
        globalparams["extendGap"] = [float(oGm.get()), float(eGp.get())]
        
        log(logs, "parameters : " + str(globalparams))
        
        SPS = 0
        bestparams = dict()
        
        for i in range(0, int(niter.get())):
            params = dict()
            for k, v in globalparams.items():
                params[k] = round(10*(v[0] + (v[1] - v[0])*random.uniform(0, 1)))/10
            
            tmpSPS = 0
            for k in SEQS.keys():
                scorer = aminoAcidScorer(str(sName.get()), params)
                tmpSPS += SPS_romainTuned(SEQSmsf[k], aligne_multiple(SEQS[k], scorer)) / len(SEQS)
            print(str(params), int(1000*tmpSPS)/1000)
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
    
    # Création d'un widget Canvas (zone graphique) pour les messages log
    Largeur = 600
    Hauteur = 320
    logger = tk.Canvas(Mafenetre, width = Largeur, height =Hauteur, bg ='white')
    logger.pack(side = tk.BOTTOM, padx = 5, pady = 5)
    
    Mafenetre.mainloop()
    
launchInterface()