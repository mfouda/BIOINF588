import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo

def Ploop1(seq):
    for i in range (0, len(seq)):
        if ((len(seq)>i+8) and (seq[i+5]=="G") and (seq[i+6]=="K") and (seq[i]=="A" or seq[i]=="G") and (seq[i+7]=="S" or seq[i+7]=="T")):
            return(True)
    return(False)

#seq=Seq("AFAKKAGAAAAGKSASSKSS")
#
#Ploop1(seq)