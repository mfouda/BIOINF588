import numpy as np

def bloc(a, b):
    return(0,a)

print(np.isnan(np.nan))

def aligne (ListOfBlocs):
    cout =0;
    list = []
    len  = len(list)
    matrix = np.zeros(len,len)
    matrix.fill(np.nan)

    #on remplit avec les distances initiales
    for i in range(len):
        for j in range(i):
            matrix[i, j] = bloc(list[i], list[j])[0]

    for iter in range(len):
        #trouve la distance minimale et les vecteurs qui la réalisent : x et y
        argmax = np.argmax(matrix)
        x, y = argmax % len, argmax//len
        min , max = min(x, y), max(x, y)

        #réalise la fusion de x et y, remplace le plus petit par la fusion et le plus grand sort du tableau (est remplacé par des nan)
        result = bloc(list[x], list[y])
        list[min] = result[1]
        matrix[max,:].fill(np.nan)
        matrix[:,max].fill(np.nan)
        for j in range(len):
            #on remplit la ligne et la colonne de min la où il n'y a pas de nan
            if not (np.isnan(matrix[min, j])):
                matrix[min,j] = bloc (list[min],list[j])[0]
            if not (np.isnan(matrix[j, min])):
                matrix[j,min] = bloc (list[j],list[min])[0]
        #actualise le cout
        cout += result[0]
    return (cout,list[0])