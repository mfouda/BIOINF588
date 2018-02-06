import numpy as np
from Bio.Seq import Seq
import Bio.SubsMat.MatrixInfo
from blocs import bloc
import random

def autoGenSeq(num, size, percentage):
    AA = ['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T']
    seqs = []
    seqs += [Seq('')]
    for i in range(0, size):
        seqs[0] += AA[random.randint(0, 19)]
    for k in range(1, num):
            seqs += [Seq('')]
            for i in range(0, size):
                if(random.uniform(0, 1) < percentage):
                    seqs[k] += seqs[0][i]
                else:
                    seqs[k] += AA[random.randint(0, 19)]
                    
    for k in range(0, num):
        seqs[k] = bloc(seqs[k], 'seq'+str(k))
        
    return seqs    

def aligne_multiple(seqs, d, e):
    score_alignement = 0
    n = len(seqs)
    
    if(type(seqs[0]) != type(bloc(Seq('')))):
        list_blocs = []
        for i in range (n) :
            list_blocs += [bloc(seqs[i], 'seq'+str(i))]
    else:
        list_blocs = seqs
    del seqs
        
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

SEQS = autoGenSeq(10, 60, 0.62)
    
#SEQS = [Seq('FEQWEKHAQYCHIVHMPDDFIVGRSNIPEKCHSLLKQCHAWIYANCIQAKGHQPLDQCETPLKPDNTWAQKYKCFVHIEKFIEFYHRSMVYHFGWCGEIC'),
#        Seq('TEQHESHAQRAWKVHMPDDFIVGTSNIPKKYQVRLKQCHAWLYENCINGYCHLPGDSFEIVLAGDNGDSQEYLPFVHIGPFIEFRTRSMVYHCGWCGEIC'),
#        Seq('FRGWILHKQYCHFVFMNNYFWVHAPRSPFKCHLLLKQKHAWPVANMIQAKAHQSLDQCFAPCKPDRTMQQKYKWFVHIEMHMEFYHESMPDHKGWFGEYC'),
#        Seq('FEQWEHHATYCHGVHMFDDDTVGRHTIWEQNESLWSQCHDWIYAKPCQPRGHQPLDQGEELLKGDNTWAQKYKVTGAAVKQSEFRHRSMLAHFKGHGINC'),
#        Seq('FEPSEKHAQYCHIWQMPDGFIYGNSNIPEQFRSWLKRHWAWITMNMITVEGVQGLYYQWTWAVDWLCWAQLYKIFVVIAPGPELDTREPVYHFMWCGEWR'),
#        Seq('KESFEGGAENDCIFCMPDDRISNHSNIPAKCHSLLMQCHHPIYDYCINRKGHQPADQYENWLETDNTWAAKYKCFVHPYKILTFWGEAMVYHFGWCVEIC'),
#        Seq('FEQWYKHAKYCHIVHMQDDFFVGRDEVPPKCHSLTLDNHAWIDPNCIQATGHKPEDRQVRPLKPDVTWAAKYKKKVYIRVFHSFTHRSMPRHFGWCGEIC'),
#        Seq('FEQWEKHAQYVHIVHMQDGFIVFGYNIPQKCQSGLRLEHSWHYFNCIQSKGHQPRQACEDPLKDDNTVAQIYRIRVYIDKFIEFNHRSMGYHFYWKGEIC'),
#        Seq('REQHHKHACYCHNEHANDDNICPHYNISECCHSELKQCHRRIYWNKKQGQMHFLLDQFETPLVPDNTTFQVYNMPVNIEGFIEFYHHHMVMVNSWCGTIC'),
#        Seq('FEQIEDFAYYAHIVHMFDDFIGGRSVICKKCHSLRRQFHAWQRANCIRAKSEYPLAQCETPLKPDNTWAWGCKYNVHIEKKFEFYHYSMVYHTIWSGMIC')]

#bloc1 = aligne_multiple([Seq("WWWAATCGTWWRTTCCTATATATCWWWTWTWWCTATCTACWWUWTWCTATZCTUAIYVWTCWA"),
#                         Seq("WWWAATCGTWWRTTCCTATUUATATCWWWTWTWWCTATCTACWWUWTUUWCTATZCTUAIYVWTCWA"),
#                         Seq("WWWAATCUGTWWRTTCCTATATATCWWWTWTWWCTUUUATCTACWWUWUUTWCTATZCTUAIYVWTCWA")],
#                        11, 1)

bloc1 = aligne_multiple(SEQS, 6, 1)
bloc1.show()