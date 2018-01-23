import numpy as np
import blocs


print(np.isnan(np.nan))

def aligne_multiple (seqs):
    score_alignement =0;
    n  = len(seqs)
    list_blocs = []
    for i in range (n) :
        new_bloc = blocs.bloc()
        new_bloc.init(seqs[i])
        list_blocs += [new_bloc]
    scores = np.zeros(n,n)
    scores.fill(np.nan)

    #on remplit avec les scores initiaux
    for i in range(n):
        for j in range(i):
            scores[i, j] = list_blocs[i].score(list_blocs[j])
    for i in range(n):
        #trouve le score maximal et les vecteurs qui le réalisent : x et y
        argmax = np.argmax(scores)
        x, y = argmax % n, argmax//n
        mini , maxi = min(x, y), max(x, y)

        #réalise la fusion de x et y, remplace le plus petit par la fusion et le plus grand sort du tableau (est remplacé par des nan)
        s = list_blocs[mini].add(list_blocs[maxi])
        scores[maxi,:].fill(np.nan)
        scores[:,maxi].fill(np.nan)
        for j in range(n):
            #on remplit la ligne et la colonne de min la où il n'y a pas de nan
            if not (np.isnan(scores[mini, j])):
                scores[mini,j] = list_blocs[mini].score(list_blocs[j])
            if not (np.isnan(scores[j, mini])):
                scores[j,mini] = list_blocs[j].score(list_blocs[mini])
        #actualise le cout
        score_alignement += s
    return (score_alignement,list_blocs[0])