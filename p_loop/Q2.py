# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:51:51 2018

@author: romai
"""

import Bio.PDB as pdb

def getInfoFromName(name):
    pdbl = pdb.PDBList()
    pdbl.retrieve_pdb_file(name)
    #http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
    
    parser = pdb.MMCIFParser()
    structure = parser.get_structure(name, 'ga/'+str(name)+'.cif')
    
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
                
getInfoFromName('2GAA')
print(22*8*2**22/1000000)