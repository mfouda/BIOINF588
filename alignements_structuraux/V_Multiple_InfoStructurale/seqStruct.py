# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:56:19 2018

@author: romai
"""

import numpy as np
from Bio.Seq import Seq
import Bio.PDB as pdb        

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
def get_num(residue):
    residue = str(residue)
    num = ""
    for i in range(0, len(residue)):
        if(RepresentsInt(residue[i])):
            num += residue[i]
    return int(num)
    
def convert_name_AA(aa):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return(d[aa])
    
def get_info(path):
    parser = pdb.PDBParser()
    structure = parser.get_structure(path, path)
    
    dico = dict()
    dico['name'] = path
    dico['num Model'] = 0
    dico['num Chain'] = 0
    dico['num Residue'] = 0
    dico['num Atom'] = 0

    baryres = [0, 0, 0]
    baryatm = [0, 0, 0]
    for model in structure:
        dico['num Model'] += 1
        for chain in model:
            dico['num Chain'] += 1
            for residue in chain:
                if(residue.get_resname() != "HOH"):
                    dico['num Residue'] += 1
                    #print("Amino acid", get_num(residue), convert_name_AA(residue.get_resname()))
                    bary_res = [0, 0, 0]
                    num_atom = 0
                    for atom in residue:
                        dico['num Atom'] += 1
                        bary_res += atom.get_coord()
                        num_atom += 1
                        baryatm += atom.get_coord()
                        #print("Atome",atom.get_coord(),atom.get_name())
                    baryres += bary_res / num_atom
                    
    dico["baryres"] = baryres / dico['num Residue']
    dico["baryatm"] = baryatm / dico['num Atom']
    return dico

def get_seq(path):
    dico = get_info(path)
    
    parser = pdb.PDBParser()
    structure = parser.get_structure(path, path)
    
    seq = dict()
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if(residue.get_resname() != "HOH"):
                    aminoacid = dict()
                    aminoacid["name"] = convert_name_AA(residue.get_resname())
                    
                    bary_res = [0, 0, 0]
                    num_atom = 0
                    for atom in residue:
                        dico['num Atom'] += 1
                        bary_res += atom.get_coord()
                        num_atom += 1
                    bary_res /= num_atom
                    #aminoacid["bary_res"] = bary_res
                    aminoacid["enfouissement"] = sum((bary_res - dico["baryres"])**2)**0.5
                    aminoacid["struct"] = "V"
                    seq[get_num(residue)] = aminoacid
    
    lines = open(path, "r").readlines()
    for line in lines:
        if(line[:6] == "HELIX "):
            start = int(line[22:25])
            end = int(line[34:37])
            for i in range(start, end+1):
                seq[i]["struct"] = "H"
                
        if(line[:6] == "SHEET "):
            start = int(line[23:26])
            end = int(line[34:37])
            for i in range(start, end+1):
                seq[i]["struct"] = "F"

    return seq


class seqStruct:
    
    def __init__(self, path):
        self.name = path[:-4]
        self.seq = []
        dico = get_seq(path)
        for cle in dico:
            self.seq += [dico[cle]]
        
    def getAminoAcid(self, i):
        return self.seq[i]
    
    def getSequence(self):
        return self.seq
    
    def getName(self):
        return self.name
    
    def printSeq(self):
        strseq = ""
        for i in range(0, len(self.getSequence())):
            strseq += self.getAminoAcid(i)["name"]
        return strseq
        
path = "2byg.pdb"
seq = seqStruct(path)
print(seq.getAminoAcid(2))
print(seq.name)