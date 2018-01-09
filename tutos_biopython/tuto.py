# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 16:10:23 2018

@author: romai
"""

import Bio.SeqIO

handle = open("ls_orchid.fasta")

i = 0
for seqrec in Bio.SeqIO.parse(handle, "fasta"):
    if(i == 0):
        print(seqrec.id)
        s = seqrec.seq
        print(s)
        print(len(s))
        
        gc = 0
        for i in range(0, len(s)):
            gc += int(s[i] == "G" or s[i] == "C")
            
        print("Taux de GC =",gc/len(s))
        
    i += 1
    
handle.close()
