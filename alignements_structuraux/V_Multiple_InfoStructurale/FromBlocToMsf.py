from blocs import*
import math

def convertBlocMsf(myBloc,name):
    n=myBloc.getNbSeqs
    l=myBloc.getLength
    nom= name+".txt"
    seqNames=myBloc.getNames()
    f = open(nom,'w')
    
    ##entête
    f.write('PileUp /n/n/n/n')
    for i in range(l):
        f.write('Name: ' + seqNames[i] + ' Len: ' + l + ' Check: ' + 4681 + ' Weight:  10.0 /n')
    f.write('/n // /n /n /n /n')
    
    ##écriture de séquence proprement dite
    nbreLigneMsf= math.ceil(l/50) 
    for i in range(nbreLigneMsf):
        for j in range(n):
            f.write(seqNames[j]+ '      ')
            k=i*50 #écrire par bloc de 50
            compteur=0 #ajouter un espace tous les dix caractères
            while ((k<(i+1)*50) and k<l):
                if myBloc.getSeq[j][k]=="-" or myBloc.getSeq[j][k]=="+": #remplacer les gaps par .
                    f.write('.')
                else:
                    f.write(myBloc.getSeq[j][k])
                compteur+=1
                if (compteur==9):
                    compteur=0
                    f.write(' ')
                k=k+1
                
    
    f.close()
    
seq = seqStruct()
l=[]
cha="WWAGCATTTGGCTGG"
for c in cha:
    dico={}
    dico["name"]=c
    l+=[dico]
seq.setSequence(l)
bloc1=bloc(seq)
#bloc1.show()

seq2 = seqStruct()
l2=[]
cha2="AGATGACTACCCT"
for c in cha2:
    dico2={}
    dico2["name"]=c
    l2+=[dico2]
seq2.setSequence(l2)
bloc2=bloc(seq2)
bloc2.show()

import score
scorer1 = score.aminoAcidScorer("blosum62", dict({"openGap" : 6, "extendGap" : 1}))
score1 = bloc1.add(bloc2,scorer1)
convertBlocMsf(bloc1,testmsf)

    

