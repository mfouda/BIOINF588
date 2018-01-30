import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo
import pickle as pk
from Bio import SeqIO



def Has_Ploop(seq):
    for i in range (0, len(seq)):
        if ((len(seq)>i+7) and (seq[i+5]=="G") and (seq[i+6]=="K") and (seq[i]=="A" or seq[i]=="G") and (seq[i+7]=="S" or seq[i+7]=="T")):
            return(True)
    return(False)

def Count_Ploop():
    dicoProt = pk.load(open( "dicoPROT.p", "rb" ))
    res = []
    for k, v in dicoProt.items():
        if Has_Ploop(v):
            res.append([k,len(v)])   
    return(res)
    

print(Count_Ploop())