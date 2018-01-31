# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:56:19 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.PDB as pdb


structure = pdb.MMCIFParser().get_structure('name', '2GAA.cif')


name = '2GAA'
dico = dict()
dico['name'] = name
dico['numM'] = 0
dico['numC'] = 0
dico['numR'] = 0
dico['numA'] = 0

for model in structure:
    dico['numM'] += 1
    for chain in model:
        dico['numC'] += 1
        for residue in chain:
            dico['numR'] += 1
            for atom in residue:
                dico['numA'] += 1
print(dico)


dssp = pdb.DSSP(structure[0], '2GAA.pdb')