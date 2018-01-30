# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:13:15 2018

@author: ASUS

ProteinAnalysis
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pickle as pk

def histo(list_of_seqs):
    dict_tot = {}
    for seqs in list_of_seqs : 
        dict_temp = ProteinAnalysis.get_amino_acids_percent(ProteinAnalysis(seqs))
        for key in dict_temp.keys():
            dict_tot[key] = dict_tot[key] + dict_temp[key]
    for key in dict_tot.keys():
            dict_tot[key] = dict_tot[key]/(float)(len(list_of_seqs))
    return dict_tot  

def histo(dico):
    return 1 #histo(list(dico.values()))
dict_temp = pk.load(open( "dicoPLOOP.p", "rb" ))
L = list(dict_temp.values())
L_def  = [str([i]) for i in L]
print(L)


dict_ploop = histo(dict_temp)
dict_tot = histo(pk.load(open( "dicoPROT.p", "rb" )))
