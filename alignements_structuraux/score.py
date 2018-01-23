import numpy as np

def aligne(a, b):
    return(0,a)

print(np.isnan(np.nan))

def aligne_multiple (blocs):
    score =0;
    n  = len(blocs)
    distances = np.zeros(n,n)
    distances.fill(np.nan)

    #on remplit avec les distances initiales
    for i in range(n):
        for j in range(i):
            distances[i, j] = aligne(blocs[i], blocs[j])[0]

    for i in range(n):
        #trouve la distance minimale et les vecteurs qui la réalisent : x et y
        argmax = np.argmax(distances)
        x, y = argmax % n, argmax//n
        mini , maxi = min(x, y), max(x, y)

        #réalise la fusion de x et y, remplace le plus petit par la fusion et le plus grand sort du tableau (est remplacé par des nan)
        result = aligne(blocs[x], blocs[y])
        blocs[mini] = result[1]
        distances[maxi,:].fill(np.nan)
        distances[:,maxi].fill(np.nan)
        for j in range(n):
            #on remplit la ligne et la colonne de min la où il n'y a pas de nan
            if not (np.isnan(distances[mini, j])):
                distances[mini,j] = aligne(blocs[mini],blocs[j])[0]
            if not (np.isnan(distances[j, mini])):
                distances[j,mini] = aligne(blocs[j],blocs[mini])[0]
        #actualise le cout
        score += result[0]
    return (score,blocs[0])