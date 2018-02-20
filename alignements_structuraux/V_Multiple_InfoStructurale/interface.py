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
        ButtonM['state'] = 'disable'
        ButtonA['state'] = 'disable'
        ButtonAps['state'] = 'disable'
        
    def Effacer():
        """ Efface la zone graphique """
        logger.delete(tk.ALL)
        
    def useTFA():
        global SEQS, params
        filename = tkfd.askopenfilename(initialdir = "../RV11/", title="Ouvrir un .tfa", filetypes=[('tfa files','.tfa'),('all files','.*')])
        print(filename)
        
        log(logs, "Recherche de séquences dans le .tfa")
        SEQS = []
        lines = open(filename, 'r').readlines()
        #lines = [lines[i][:-2] for i in range(0, len(lines))]
        for line in lines:
            if(line[0] == ">"):
                print(line[1:5])
                if(line[1:5] + ".pdb" not in os.listdir("pdb/")):
                    log(logs, "Downloading " + line[1:5] + ".pdb ...")
                    path = "https://files.rcsb.org/download/" + line[1:5] + ".pdb"
                    urllib.urlretrieve(path, "pdb/" + line[1:5] + ".pdb")
                else:
                    log(logs, "Protein " + line[1:5] + ".pdb already downloaded before")
                SEQS += [seqStruct("pdb/" + line[1:5] + ".pdb")]
        #.split("<")
        #print(lines)
#        for i in range(0, int(nM.get())):
#            SEQS += [SEQS[0].mutate(float(Mtx.get()))]
            
        params = dict()
        ButtonA['state'] = 'normal'
        ButtonAps['state'] = 'normal'
        showSEQS(SEQS)
    
    def log(logs, s):
        Effacer()
        logs += ["[" + str(datetime.datetime.now())[:-7] + "]    " + s]
        nrows = int(Hauteur / 10)
        for i in range(0, min(nrows, len(logs))):
            logger.create_text(5, 5 + 1.1*(i / nrows) * (Hauteur - 10), text = logs[-i-1], anchor = tk.NW, font=("Helvetica", 9))
            
    def showP():
        global SEQS
        getP()
        SEQS = [seqStruct("pdb/" + strPName.get() + ".pdb")]
        ButtonM['state'] = 'normal'
        log(logs, "Protéine "+str(strPName.get())+" initialisée")
        showProtein(strPName.get())
        
    def getP():
        if(strPName.get() + ".pdb" not in os.listdir("pdb/")):
            log(logs, "Downloading " + strPName.get() + ".pdb ...")
            path = "https://files.rcsb.org/download/" + strPName.get() + ".pdb"
            urllib.urlretrieve(path, "pdb/" + strPName.get() + ".pdb")
        else:
            log(logs, "Protein " + strPName.get() + ".pdb already downloaded before")
            
    def mutate():
        global SEQS, params
        log(logs, "Mutation "+str(nM.get())+" fois avec "+Mtx.get()+"% de similitude")
        for i in range(0, int(nM.get())):
            SEQS += [SEQS[0].mutate(float(Mtx.get()))]
        ButtonA['state'] = 'normal'
        ButtonAps['state'] = 'normal'
        params = dict()
        showSEQS(SEQS)
    
    def addParams():
        global params
        params[str(key.get())] = float(value.get())
        log(logs, "Adding entry in params  {" + str(key.get()) + " : " + str(value.get())+"}")
        
    def alignate():
        global SEQS
        global params
        params["openGap"] = float(oG.get())
        params["extendGap"] = float(eG.get())
        scorer = aminoAcidScorer(str(sName.get()), params)
        log(logs, "scorer : " + scorer.getName() + " " + str(scorer.getParams()))
        log(logs, "Alignement des séquences en cours ...")
        bloc = aligne_multiple(SEQS, scorer)
        showBloc(bloc)
            
    # Création de la fenêtre principale (main window)
    Mafenetre = tk.Tk()
    Mafenetre.title('Protein alignment')
    
    #Widget Protein
    FramePName = tk.Frame(Mafenetre,borderwidth=2,relief=tk.GROOVE)
    FramePName.pack(side=tk.LEFT,padx=2,pady=2)
    LabelP = tk.Label(FramePName,text="Proteine à analyser")
    LabelP.pack(padx=2,pady=2)
    strPName = tk.StringVar()
    strPName.set("2byg")
    TextVarP = tk.Entry(FramePName, textvariable=strPName, width=15)
    TextVarP.pack(padx=2,pady=2)
    ButtonP = tk.Button(FramePName,text="Afficher la proteine",fg='navy',command=showP)
    ButtonP.pack(padx=2,pady=2)
    
    ButtonTFA = tk.Button(FramePName,text="Utiliser .tfa",fg='navy',command=useTFA)
    ButtonTFA.pack(padx=2,pady=2)
    
    
    #Widget Mutation
    FrameM = tk.Frame(Mafenetre,borderwidth=2,relief=tk.GROOVE)
    FrameM.pack(side=tk.LEFT,padx=2,pady=2)
    LabelM = tk.Label(FrameM,text="Taux de mutation et nombre de mutants")
    LabelM.pack(padx=2,pady=2)
    Mtx = tk.StringVar()
    Mtx.set(0.62)
    TextVarM = tk.Entry(FrameM, textvariable=Mtx, width=10)
    TextVarM.pack(padx=2,pady=2)
    nM = tk.StringVar()
    nM.set(2)
    TextVarnM = tk.Entry(FrameM, textvariable=nM, width=10)
    TextVarnM.pack(padx=2,pady=2)
    ButtonM = tk.Button(FrameM,text="Muter",fg='navy',command=mutate, state=tk.DISABLED)
    ButtonM.pack(padx=2,pady=2)
    
    
    #Widget Scorer
    FrameS = tk.Frame(Mafenetre,borderwidth=2,relief=tk.GROOVE)
    FrameS.pack(side=tk.LEFT,padx=2,pady=2)
    LabelS= tk.Label(FrameS,text="Scorer")
    LabelS.pack(padx=2,pady=2)
    sName = tk.StringVar()
    sName.set("blosum62")
    TextVarsName = tk.Entry(FrameS, textvariable=sName, width=10)
    TextVarsName.pack(padx=2,pady=2)
    oG = tk.StringVar()
    oG.set(6)
    TextVaroG = tk.Entry(FrameS, textvariable=oG, width=10)
    TextVaroG.pack(padx=2,pady=2)
    eG = tk.StringVar()
    eG.set(1)
    TextVareG = tk.Entry(FrameS, textvariable=eG, width=10)
    TextVareG.pack(padx=2,pady=2)
    
    FrameAps = tk.Frame(FrameS,borderwidth=2,relief=tk.GROOVE)
    FrameAps.pack(padx=2,pady=2)
    key = tk.StringVar()
    key.set("key")
    TextVarkey = tk.Entry(FrameAps, textvariable=key, width=10)
    TextVarkey.pack(padx=2,pady=2)
    value = tk.StringVar()
    value.set(1)
    TextVarvalue = tk.Entry(FrameAps, textvariable=value, width=10)
    TextVarvalue.pack(padx=2,pady=2)
    ButtonAps = tk.Button(FrameAps,text="Add parameter",fg='navy',command=addParams, state=tk.DISABLED)
    ButtonAps.pack(padx=2,pady=2)
    
    ButtonA = tk.Button(FrameS,text="Aligner",fg='navy',command=alignate, state=tk.DISABLED)
    ButtonA.pack(padx=2,pady=2)
    
    
    # Création d'un widget Button (bouton Effacer)
    BoutonReset = tk.Button(Mafenetre, text ='Reset', command = Reset)
    BoutonReset.pack(side = tk.LEFT, padx = 5, pady = 5)
    
    # Création d'un widget Button (bouton Quitter)
    BoutonQuitter = tk.Button(Mafenetre, text ='Quitter', command = Mafenetre.destroy)
    BoutonQuitter.pack(side = tk.LEFT, padx = 5, pady = 5)
    
    # Création d'un widget Canvas (zone graphique) pour les messages log
    Largeur = 600
    Hauteur = 320
    logger = tk.Canvas(Mafenetre, width = Largeur, height =Hauteur, bg ='white')
    logger.pack(padx = 5, pady = 5)
    
    Mafenetre.mainloop()
    
launchInterface()