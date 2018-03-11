# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 13:59:18 2018

@author: Jasmine
"""

from blocalignments import aligne_multiple
from blocs import bloc
from seqStruct import seqStruct
from Bio import SeqIO


def parse_tfa (filename) :
    seqs = []
    lines = open(filename, 'r').readlines()
    n = len(lines)
    i = 0
    while i<n :
        seq = ''
        i += 2
        while i<n and not lines[i] == '\n' :
            l = len(lines[i])
            seq += lines[i][:l-1]
            i += 1
        seqs.append(seq)
        while i<n and not lines[i][0] == '>' :
            i+=1
    return seqs

#transforme un fichier msf (alignement de référence) en dictionnaire de SeqStruct
def msftoDict (filename) :
    seqs = []
    lines = open(filename, 'r').readlines()
    n = len(lines)
    nb_seqs = 0
    i = 6
    while lines[i][:5] == ' Name' :
        nb_seqs += 1
        i += 1
    for k in range (nb_seqs) :
        seqs.append('')
    i += 5
    res = dict()
    last = dict()
    
    while (i + nb_seqs - 1) < n :
        for k in range (nb_seqs) :
            words = lines[i+k].split()
            name = words[0]
            if name not in res :
                res[name] = seqStruct()
                res[name].setName(name)
                last[name] = 0
                
            for word in words[1:] :
                for aa in word :
                    d = dict()
                    d["name"] = aa
                    if(aa != "."):
                        d["id"]=last[name]
                        last[name] += 1        
                    res[name].addAminoAcidAfter(d)
        i += nb_seqs + 2

    return res


#st est un string représentant une séquence sans gap
def strToSeqStruct(name, st) :
    seq = seqStruct()
    seq.setName(name)
    i = 0
    for c in st :
        d = dict()
        d["name"] = c
        d["id"] = i
        seq.addAminoAcidAfter(d)
        i += 1
    return seq


def TC (ref_msf,our_result):
    seqs = [s for s in our_result.getSeqs()]
    return score_TC_computer(seqs, ref_msf)

def TC_romainTuned (ref_msf,our_result):
    return score_TC_computer(our_result, ref_msf)


#seqs stocke notre alignement sous la forme d'une liste de seqStruct
#dico stocke l'alignement de référence sous la forme d'un dico de seqStruct
# /!\ c'est le contraire pour SPS
def score_TC_computer (seqs, dico) :
    nb_seq = len(seqs)
    nb_col = seqs[0].getLength()
    nb_col_ref = dico[seqs[0].getName()].getLength()
    i2 = [0]*nb_seq #i2[k] contient la position du premier aa non regardé de la séquence k de référence
    score = 0

    for i in range(nb_col) :
        res = True
        j = 0
        
        while j< nb_seq and 'id' not in seqs[j].getAminoAcid(i):
            j += 1
        seq_ref = dico[seqs[j].getName()]
        while i2[j] < nb_col_ref and 'id' not in seq_ref.getAminoAcid(i2[j]) :
            i2[j] += 1
        #j est le numero de la première séquence non vide de la colonne i (=n'ayant pas de gap en position i) dans notre alignement
        #i2[j] est le numero de la colonne de référence correspondante
        
        #dans cette boucle on regarde si les aa de la colonne i de notre alignement sont tous situés dans la colonne i2[j] de l'alignement de référence
        for k in range (j+1,nb_seq) :
            our_aa = seqs[k].getAminoAcid(i)
            if 'id' in our_aa :
                aa_ref = dico[seqs[k].getName()].getAminoAcid(i2[j])
                while i2[k]< nb_col_ref and 'id' not in dico[seqs[k].getName()].getAminoAcid(i2[k]) :
                    i2[k] += 1
                i2[k] += 1
                res = res*('id' in aa_ref and aa_ref['id']==our_aa['id'])
        i2[j] += 1
        
        score += int(res)
    
    return score/nb_col


def test_TC(filename, scorer) :
    fasta_seqs = SeqIO.parse(open('../RV11/' + filename + '.tfa'), 'fasta')
    seqs = []
    for fs in fasta_seqs :
        seqs.append(strToSeqStruct(fs.id, str(fs.seq)))
    our_bloc = aligne_multiple(seqs, scorer)
    ref = msftoDict('../RV11/' + filename + '.msf')
    print(TC(ref, our_bloc))

#import score
#scorer1 = score.aminoAcidScorer("blosum62", dict({"openGap" : 11, "extendGap" : 1}))
##scorer2 = aminoAcidScorer("blosum62mixte", dict({"openGap" : 6, "extendGap" : 1}))
#test_TC('BB11001', scorer1)

