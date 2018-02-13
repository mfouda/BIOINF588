# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 14:47:55 2018
@author: romai
"""

##################
##### p-loop #####
##################

# 1/ Détection de signatures sur un génome complet
# 1.1/ Identification du motif

# Question 1 : Fonction qui télécharge un fichier

import urllib.request as urllib

def downloadPrositeData(path):
    urllib.urlretrieve(path, "data/prosite.dat")
    
########## Executer Q1 ###########
#path = "ftp://ftp.expasy.org/databases/prosite/prosite.dat"
#downloadPrositeData(path)
##################################

# Question 2 : Recherche du motif du p-loop dans le fichier

from itertools import product

def findSomethingInFile(path, smthg):
    lines = open(path, "r").readlines()
    
    L = [(c, c.upper()) if not c.isdigit() else (c,) for c in smthg.lower()]
    L = ["".join(item) for item in product(*L)]
    
    for i in range(0, len(lines)):
        find = False
        for l in L:
            if(lines[i].find(l) != -1):
                find = True
        if(find):
            print(lines[i], lines[i+1])

########## Executer Q2 ###########
#path = "data/prosite.dat"
#findSomethingInFile(path, "p-loop")
##################################

# 1.2/ Lecture du format FASTA

# Question 1 : Fonction qui télécharge un fichier fasta

import gzip as gzip

def downloadFasta(path, dest):
    urllib.urlretrieve(path, dest+".gz")
    f = gzip.open(dest+".gz", 'rb')
    pk.dump(f.read(), open( dest, "wb" ) )
    
########## Executer Q1 ###########
#path = "https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz"
#downloadFasta(path, "data/orf_coding_all.fasta")
##################################

# Question 2 : Lire FASTA puis stocker au format approprié

from Bio import SeqIO
import pickle as pk
import numpy as np

def fromFastaARNToPickleARN(path):
    fasta_sequences = SeqIO.parse(open(path),'fasta')

    dicoARN = dict()
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        dicoARN[name] = sequence
    
    return dicoARN

########## Executer Q2 ###########
#path = "data/orf_coding_all.fasta"
#dicoARN = fromFastaARNToPickleARN(path)
#pk.dump(dicoARN, open( "pickleObjects/dicoARN.p", "wb" ) )
#print(dicoARN)
##################################

# Question 3 : Traduire ARN en acide aminé

from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

def fromPickleARNtoPickleProtein(path):
    dicoARN = pk.load(open( path, "rb" ))
    
    dicoPROT = dict()
    for k, v in dicoARN.items():
        dicoPROT[k] = Seq(v, generic_rna).translate()
        
    return dicoPROT
    
########## Executer Q3 ###########
#path = "pickleObjects/dicoARN.p"
#dicoPROT = fromPickleARNtoPickleProtein(path)
#pk.dump(dicoPROT, open( "pickleObjects/dicoPROT.p", "wb" ) )
#print(dicoPROT)
##################################

# Question 4 : Comparer le fichier obtenu avec un autre fichier

def comparePickleWithFasta(picklePath, fastaPath):
    print("Comparing the two files ...")
    fasta_sequences = SeqIO.parse(open(fastaPath),'fasta')
    dicoPROT = pk.load(open( picklePath, "rb" ))
    
    areEquals = True
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if(dicoPROT[name] != sequence):
            areEquals = False
            print(dicoPROT[name], sequence)
            print("")
            
    if(areEquals):
        print("Les deux fichiers sont égaux")
    else:
        print("Les deux fichiers sont différents")
    


########## Executer Q3 ###########
path = "https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz"
downloadFasta(path, "data/orf_trans_all.fasta")
comparePickleWithFasta("pickleObjects/dicoPROT.p", "data/orf_trans_all.fasta")
##################################
    
# 1.3/ Identification des protéines à p-loop

# Question 1 : Recherche des p-loops
    
def hasPloop(seq):
    for i in range (0, len(seq)):
        if ((len(seq)>i+7)
        and (seq[i+5]=="G")
        and (seq[i+6]=="K")
        and (seq[i]=="A" or seq[i]=="G")
        and (seq[i+7]=="S" or seq[i+7]=="T")):
            return(True)
    return(False)

def findPloop(path):
    dicoProt = pk.load(open(path, "rb" ))
    
    dicoPLOOP_PROT = dict()
    for k, v in dicoProt.items():
        if hasPloop(v):
            dicoPLOOP_PROT[k] = v
            
    return dicoPLOOP_PROT

########## Executer Q1 ###########
#path = "pickleObjects/dicoPROT.p"
#dicoPLOOP_PROT = findPloop(path)
#pk.dump(dicoPLOOP_PROT, open( "pickleObjects/dicoPLOOP_PROT.p", "wb" ) )
#print(len(dicoPLOOP_PROT))
##################################
    
# Question 2 : Histogramme de la composition en acide aminé des proteines à p-loop
    
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def aminoAcidComposition(path):
    dicoProt = pk.load(open(path, "rb" ))
    
    hist = dict()
    for k, v in dicoProt.items():
        hist_temp = ProteinAnalysis.get_amino_acids_percent(ProteinAnalysis(str(v)))
        for key in hist_temp.keys():
            if(key in hist):
                hist[key] = hist[key] + 100 * hist_temp[key] / (float)(len(dicoProt))
            else :
                hist[key] = 100 * hist_temp[key] / (float)(len(dicoProt))
    
    return hist

def histogramFromCompos(comps, labels):
    plt.figure(figsize = (16, 8))
    
    delta = 3 * (np.linspace(0, 1, len(comps)) - [0.5]*len(comps)) / 4
    
    ymax = 0
    for i in range(0, len(comps)):
        comp = comps[i]
        label = labels[i]
        if(len(comps) != 1):
            X = np.arange(len(comp)) + [delta[i] - i*0.75 / (2*(len(comps) - 1))]*20
        else:
            X = np.arange(len(comp))
        plt.bar(X, comp.values(), align='edge', width = 0.75 / len(comps), label = label)
        ymax = max(ymax, max(comp.values())*1.05)
        
    plt.ylim(0, ymax)  
    plt.xticks(np.arange(20), comp.keys())
    plt.title("Amino acid composition")
    plt.xlabel("Amino acid")
    plt.ylabel("Proportion")
    plt.legend(loc = "best")
    plt.show()

########## Executer Q2 ###########
#path = "pickleObjects/dicoPLOOP_PROT.p"
#histPLOOP_PROT = aminoAcidComposition(path)
#path = "pickleObjects/dicoPROT.p"
#histPROT = aminoAcidComposition(path)
#histogramFromCompos([histPROT, histPLOOP_PROT], ["prot","ploop_prot"])
##################################

# 1.4/ Enzymes de restriction

# Question 1 : Donner sites de restriction de certaines enzymes

from Bio import Restriction

########## Executer Q1 ###########
#print("Site de restriction de l'enzyme", Restriction.EcoRI, ":", Restriction.EcoRI.site)
#print("Site de restriction de l'enzyme", Restriction.XhoI, ":", Restriction.XhoI.site)
#print("Site de restriction de l'enzyme", Restriction.TaqI, ":", Restriction.TaqI.site)
##################################

# Question 2 : Trouver les sites de restriction parmi les proteines à p-loop

def findRestrictionSites(pathARN, pathPLOOP_PROT):
    dicoPLOOP_PROT = pk.load(open( pathPLOOP_PROT, "rb" ))
    dicoARN = pk.load(open( pathARN, "rb" ))

    dicoEcoRI = dict()
    dicoXhoI = dict()
    dicoTaqI = dict()

    for k, v in dicoPLOOP_PROT.items():
        dicoEcoRI[k] = Restriction.EcoRI.search(Seq(dicoARN[k]))
        dicoXhoI[k] = Restriction.XhoI.search(Seq(dicoARN[k]))
        dicoTaqI[k] = Restriction.TaqI.search(Seq(dicoARN[k]))
    
    return dicoEcoRI, dicoXhoI, dicoTaqI

########## Executer Q2 ###########
#pathARN, pathPLOOP_PROT = "pickleObjects/dicoARN.p", "pickleObjects/dicoPLOOP_PROT.p"
#dicoEcoRI, dicoXhoI, dicoTaqI = findRestrictionSites(pathARN, pathPLOOP_PROT)
#pk.dump(dicoEcoRI, open( "pickleObjects/dicoEcoRI.p", "wb" ) )
#pk.dump(dicoXhoI, open( "pickleObjects/dicoXhoI.p", "wb" ) )
#pk.dump(dicoTaqI, open( "pickleObjects/dicoTaqI.p", "wb" ) )
##################################
    
# 2/ Analyse de séquence et structure 3D
    
# Question 1 : Télécharger fichier PDB a partir du code

def downloadPDBFromName(name):
    path = "https://files.rcsb.org/download/" + name + ".pdb"
    urllib.urlretrieve(path, "pdbFiles/" + name + ".pdb")
    
########## Executer Q1 ###########
#name = "2GAA"
#downloadPDBFromName(name)
##################################
    
# Question 2 : Utiliser le parser BioPython pour récuperer des infos nécessaires
    
import Bio.PDB as pdb

def getInfoFromPDBFile(path):
    parser = pdb.PDBParser()
    structure = parser.get_structure(path, path)
    
    dico = dict()
    dico['name'] = path[-8:-4]
    dico['num Model'] = 0
    dico['num Chain'] = 0
    dico['num Residue'] = 0
    dico['num Atom'] = 0

    for model in structure:
        dico['num Model'] += 1
        for chain in model:
            dico['num Chain'] += 1
            for residue in chain:
                dico['num Residue'] += 1
                for atom in residue:
                    dico['num Atom'] += 1
                    
    return dico

########## Executer Q2 ###########
#path = "pdbFiles/2GAA.pdb"
#print(getInfoFromPDBFile(path))
##################################
    
# Question 3 :  ...
    
# ...







