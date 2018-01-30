# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:13:15 2018

@author: ASUS
"""

from Bio import Bio.ProteinAnalysis

def histo(list):
    dict_tot = {}
    for seqs in list : 
        dict = ProteinAnalysis(seqs).get_amino_acids_percent()
        for key in dict.keys():
            dict_tot[key] = dict_tot[key] + dict[key]
    for key in dict_tot.keys():
            dict_tot[key] = dict_tot[key]/(float)list.len()
    return dict_tot  