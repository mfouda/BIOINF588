import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo
from blocs import bloc  
from seqStruct import seqStruct
from score import aminoAcidScorer

def aligne_multiple(seqs, scorer):
    score_alignement = 0
    n = len(seqs)
    
    list_blocs = []
    for i in range(0, n) :
        list_blocs += [bloc(seqs[i])]
    del seqs
    
    scores = np.zeros((n, n))
    scores.fill(np.nan)
    #on remplit avec les scores initiaux
    for i in range(n):
        for j in range(i):
            scores[i, j] = list_blocs[i].alignementscore(list_blocs[j], scorer)
            
    for i in range(0, n-1):
        #trouve le score maximal et les vecteurs qui le réalisent : x et y
        argmax = np.nanargmax(scores)
        x, y = argmax % n, argmax//n
        mini , maxi = min(x, y), max(x, y)
        #réalise la fusion de x et y, remplace le plus petit par la fusion et le plus grand sort du tableau (est remplacé par des nan)
        s = list_blocs[mini].add(list_blocs[maxi], scorer)
        scores[maxi,:].fill(np.nan)
        scores[:,maxi].fill(np.nan)
        for j in range(n):
            #on remplit la ligne et la colonne de min la où il n'y a pas de nan
            if not (np.isnan(scores[mini, j])):
                scores[mini,j] = list_blocs[mini].alignementscore(list_blocs[j], scorer)
            if not (np.isnan(scores[j, mini])):
                scores[j,mini] = list_blocs[j].alignementscore(list_blocs[mini], scorer)
        #actualise le cout
        score_alignement += s
    return list_blocs[0]
    

path = "pdb/2byg.pdb"
seq = seqStruct(path)
SEQS = [seq]
for i in range(0, 3):
    SEQS += [SEQS[0].mutate(0.8)]
    
print("")
scorer1 = aminoAcidScorer("blosum62", dict({"openGap" : 4, "extendGap" : 1}))
scorer2 = aminoAcidScorer("blosum62mixte", dict({"openGap" : 6, "extendGap" : 1}))
bloc1 = aligne_multiple(SEQS, scorer1)
bloc1.show()