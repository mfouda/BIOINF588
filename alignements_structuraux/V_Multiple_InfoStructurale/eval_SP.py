from blocs import bloc

def SP(myBloc): ##https://www.cs.princeton.edu/~mona/Lecture/msa1.pdf doc pour créer Sum of Pairs
    n=myBloc.getNbSeqs
    l=myBloc.getLength
    score=0
    for i in range(l):
        col=myBloc.getCol(i)
        for j in range(n):
            for k in range(j):
                if not (col[j]==col[k] and (col[j]=="-" or col[j]=="+")):              
                    if col[j]==col[k]:
                        score +=1          
                    elif col[j]=="-" or col[j]=="+":
                        score-=2
                    else:
                        score -=1
    return(score)
