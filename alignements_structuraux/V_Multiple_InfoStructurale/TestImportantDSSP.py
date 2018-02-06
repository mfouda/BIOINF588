# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:56:19 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.PDB as pdb

name = "2byg"
parser = pdb.PDBParser()
structure = parser.get_structure(name, name + '.pdb')
repr(structure)

dico = dict()
dico['name'] = name
dico['num Model'] = 0
dico['num Chain'] = 0
dico['num Residue'] = 0
dico['num Atom'] = 0

i = 0

for model in structure:
    dico['num Model'] += 1
    for chain in model:
        dico['num Chain'] += 1
        for residue in chain:
            dico['num Residue'] += 1
            i += 1
            print(i, residue.get_resname())
            for atom in residue:
                dico['num Atom'] += 1

print(dico)