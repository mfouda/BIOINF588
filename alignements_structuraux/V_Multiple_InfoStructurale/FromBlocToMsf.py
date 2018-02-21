from blocs import*

def convertBlocMsf(myBloc,name):
    n=myBloc.getNbSeqs
    l=myBloc.getLength
    nom= name+".txt"
    seqNames=myBloc.getNames()
    f = open(nom,'w')
    f.write('PileUp /n/n/n/n')
    for i in range(l):
        f.write('Name: ' + seqNames[i] + ' Len: ' + l + ' Check: ' + 4681 + ' Weight:  10.0 /n')
    f.write('/n // /n /n /n /n')
    nbreLigneMsf= maths.roof(l/50)
    for i in range(nbreLigneMsf):
        for j in range(n):
            f.write(seqNames[j]+ '      ')
            k=i*50
            while ((k<(i+1)*50) and k<l):
                if myBloc.getSeq[j][k]=="-" or myBloc.getSeq[j][k]=="+": ##tout les 10 espace
                    f.write('.')
                else:
                    f.write()
                k=k+1
                
    
    f.close()
    

'