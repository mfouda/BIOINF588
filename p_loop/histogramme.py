# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:13:15 2018

@author: ASUS

ProteinAnalysis
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pickle as pk

def histo_of_list(list_of_seqs):
    dict_tot = {}
    for seqs in list_of_seqs : 
        dict_temp = ProteinAnalysis.get_amino_acids_percent(ProteinAnalysis(seqs))
        for key in dict_temp.keys():
            if(key in dict_tot):
                dict_tot[key] = dict_tot[key] + dict_temp[key]
            else : 
                dict_tot[key] = dict_temp[key]
            
    for key in dict_tot.keys():
            dict_tot[key] = dict_tot[key]/(float)(len(list_of_seqs))
    return dict_tot  

def histo_of_dico(dico):
    return histo_of_list(list(dico.values()))

def print_histo(list_of_histo, nombre_de_cs = 3):
    for key in list_of_histo[0].keys():
        print (key)
        for histo in list_of_histo : 
            if(key in histo):
                print (round(histo[key]*100,nombre_de_cs),'%')
            else : 
                print (0)

#prétraiter les données car elles sortent sous forme de sequences et qu'il faut des Strings
dict_of_seq_loop = pk.load(open( "dicoPLOOP.p", "rb" ))
L_of_seq_loop = list(dict_of_seq_loop.values())
L_of_string_loop  = [str(i) for i in L_of_seq_loop]

dict_of_seq_tot = pk.load(open( "dicoPROT.p", "rb" ))
L_of_seq_tot = list(dict_of_seq_tot.values())
L_of_string_tot  = [str(i) for i in L_of_seq_tot]

#on réalise les deux histogrammes
histo_ploop = histo_of_list(L_of_string_loop)
histo_tot = histo_of_list(L_of_string_tot)

#print
list_histo = [histo_ploop, histo_tot]
print_histo(list_histo)
