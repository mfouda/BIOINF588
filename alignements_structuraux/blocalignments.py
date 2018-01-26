import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo
from blocs import bloc

def aligne_multiple(seqs, d, e):
    score_alignement = 0
    n = len(seqs)
    list_blocs = []
    for i in range (n) :
        list_blocs += [bloc(seqs[i])]
    scores = np.zeros((n, n))
    scores.fill(np.nan)
    #on remplit avec les scores initiaux
    for i in range(n):
        for j in range(i):
            scores[i, j] = list_blocs[i].alignementscore(list_blocs[j], d, e)
            
    for i in range(n-1):
        #trouve le score maximal et les vecteurs qui le réalisent : x et y
        argmax = np.nanargmax(scores)
        x, y = argmax % n, argmax//n
        mini , maxi = min(x, y), max(x, y)
        
        #réalise la fusion de x et y, remplace le plus petit par la fusion et le plus grand sort du tableau (est remplacé par des nan)
        s = list_blocs[mini].add(list_blocs[maxi], d, e)
        scores[maxi,:].fill(np.nan)
        scores[:,maxi].fill(np.nan)
        for j in range(n):
            #on remplit la ligne et la colonne de min la où il n'y a pas de nan
            if not (np.isnan(scores[mini, j])):
                scores[mini,j] = list_blocs[mini].alignementscore(list_blocs[j], d, e)
            if not (np.isnan(scores[j, mini])):
                scores[j,mini] = list_blocs[j].alignementscore(list_blocs[mini], d, e)
        #actualise le cout
        score_alignement += s
        
    return list_blocs[0]

bloc = aligne_multiple([Seq("WWAACATTTGGCTGG"), Seq("WWAGCATTTAGCTGG"), Seq("WWAGCAATTGGCTGG"), Seq("AGATGACTACCCT")], 11, 1)
bloc.show()