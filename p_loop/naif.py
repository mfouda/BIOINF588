# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:55:57 2018

@author: Jasmine
"""

def pLoop_naif (seq) :
    result = False
    n = len(seq)
    for i in range(n-8) :
        if seq[i] == 'A' or seq[i] == 'G' :
            if seq[i+5] == 'G' and seq[i+6] == 'K' :
                if seq[i+7] == 'S' or seq[i+7] == 'T' :
                    result = True
    return result