# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:21:02 2018

@author: romai
"""
from Bio.Seq import Seq
from Bio import Restriction
import pickle as pk

print(Restriction.EcoRI, Restriction.EcoRI.site)
print(Restriction.XhoI, Restriction.XhoI.site)
print(Restriction.TaqI, Restriction.TaqI.site)

dicoPLOOP = pk.load(open( "dicoPLOOP.p", "rb" ))
dicoARN = pk.load(open( "dicoARN.p", "rb" ))

dicoEcoRI = dict()
dicoXhoI = dict()
dicoTaqI = dict()

for k, v in dicoPLOOP.items():
    dicoEcoRI[k] = Restriction.EcoRI.search(Seq(dicoARN[k]))
    dicoXhoI[k] = Restriction.XhoI.search(Seq(dicoARN[k]))
    dicoTaqI[k] = Restriction.TaqI.search(Seq(dicoARN[k]))
    
pk.dump(dicoEcoRI, open( "dicoEcoRI.p", "wb" ) )
pk.dump(dicoXhoI, open( "dicoXhoI.p", "wb" ) )
pk.dump(dicoTaqI, open( "dicoTaqI.p", "wb" ) )