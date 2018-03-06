from blocalignments import aligne_multiple
from blocs import bloc
from seqStruct import seqStruct
from Bio import SeqIO

def ordonne_col(j, n, aligne_ref, bloc_result) :
    col = bloc_result.getCol(j)
    col_ref = []
    col_nous = ['']*n
    for i in range (n) :
        col_ref.append(aligne_ref[i][j])
        name_seq = bloc_result.getName(i)
        k = int(name_seq[3:])
        col_nous[k] = col[i]
    return col_nous, col_ref

def parse_tfa (filename) :
    seqs = []
    lines = open(filename, 'r').readlines()
    n = len(lines)
    i = 0
    while i<n :
        seq = ''
        i += 2
        while i<n and not lines[i] == '\n' :
            l = len(lines[i])
            seq += lines[i][:l-1]
            i += 1
        seqs.append(seq)
        while i<n and not lines[i][0] == '>' :
            i+=1
    return seqs

def msftoBloc (filename) :
    seqs = []
    lines = open(filename, 'r').readlines()
    n = len(lines)
    nb_seqs = 0
    i = 6
    while lines[i][:5] == ' Name' :
        nb_seqs += 1
        i += 1
    for k in range (nb_seqs) :
        seqs.append('')
    i += 5
    res = dict()
    last = dict()
    
    while (i + nb_seqs - 1) < n :
        for k in range (nb_seqs) :
            words = lines[i+k].split()
            name = words[0]
            if name not in res :
                res[name] = seqStruct()
                res[name].setName(name)
                last[name] = 0
                
            for word in words[1:] :
                for aa in word :
                    d = dict()
                    d["name"] = aa
                    if(aa != "."):
                        d["id"]=last[name]
                        last[name] += 1        
                    res[name].addAminoAcidAfter(d)
        i += nb_seqs + 2

    return res


#st est un string représentant une séquence sans gap
def strToSeqStruct(name, st) :
    seq = seqStruct()
    seq.setName(name)
    i = 0
    for c in st :
        d = dict()
        d["name"] = c
        d["id"] = i
        seq.addAminoAcidAfter(d)
        i += 1
    return seq



#def score_C2 (filename, d, e) :
#    seqs_input = parse_tfa((filename + '.tfa'))
#    bloc_result = aligne_multiple(seqs_input, d, e)
#    aligne_ref = parse_msf(filename + '.msf')
#    
#    n = bloc_result.getNbSeqs()
#    m = len(bloc_result.getSeqs()[0])
#    score_C = 0
#    score_C2 = 0
#    
#    fichier = open(filename + '_score_C2.txt', 'w')
#    carac_per_line = 70
#    lines = ['']*(2 + m//carac_per_line)*(2*n + 8)
#    line = 0
#    carac = 0
#    
#    for j in range (m) :
#        col_nous, col_ref = ordonne_col(j, n, aligne_ref, bloc_result)
#        score_col = 0
#        
#        if carac == carac_per_line :
#            line += 2*n + 6
#            carac = 0
#        
#        for i in range (n) :
#            lines[line + i] += col_nous[i] + ' '
#            lines[line + i + n + 1] += col_ref[i] + ' '
#            score_col += (col_nous[i] == col_ref[i])
#        
#        lines[line + 2*n + 2] += str(score_col) + ' '
#        carac += 1           
#        score_C2 += score_col
#        score_C += (score_col == n)
#    
#    lines[line + 2*n + 5] += 'score_C2 = ' + str(score_C2) + ' / ' + str(m*n) + ' = ' + str(score_C2/(m*n))
#    lines[line + 2*n + 4] += 'score_C = ' + str(score_C) + ' / ' + str(m) + ' = ' + str(score_C/m)
#    
#    for x in lines :
#        fichier.write(x + '\n')
#    fichier.close()
#    
#    return score_C/m, score_C2/(m*n)

def SPS(ref_msf,our_result):
    seqs = [v for v in ref_msf.values()]
    dico = dict()
    for l in our_result.getSeqs() :
        dico[l.getName()] = l
    return score_SPS_computer(seqs, dico)

def SPS_romainTuned(ref_msf,our_result):
    dico = dict()
    for l in our_result.getSeqs() :
        dico[l.getName()] = l
    return score_SPS_computer(ref_msf, dico)

def score_SPS_computer (seqs, dico) :
    score = 0
    score_ref = 0
    nb_seq = len(seqs)
    nb_col = seqs[0].getLength()
    
    for i in range(nb_seq-1) :
        seq1 = seqs[i]
        seq2 = dico[seq1.getName()]
        j2 = 0
        for j in range(nb_col) :
            if 'id' in seq1.getAminoAcid(j):
                while 'id' not in seq2.getAminoAcid(j2) :
                    j2 += 1
                col = []
                for k in range(i+1, nb_seq):
                    col += [(seqs[k].getAminoAcid(j), seqs[k].getName())]
                for aa in col :
                    if 'id' in aa[0] :
                        score_ref += 1
                        aa2 = dico[aa[1]].getAminoAcid(j2)
                        if 'id' in aa2 :
                            score += (aa[0]["id"] == aa2["id"])          
        
    return score / score_ref

def test_SPS(filename, scorer) :
    fasta_seqs = SeqIO.parse(open('../RV11/' + filename + '.tfa'), 'fasta')
    seqs = []
    for fs in fasta_seqs :
        seqs.append(strToSeqStruct(fs.id, str(fs.seq)))
    our_bloc = aligne_multiple(seqs, scorer)
    print(SPS(msftoBloc('../RV11/' + filename + '.msf'), our_bloc))

#import score
#scorer1 = score.aminoAcidScorer("blosum62", dict({"openGap" : 6, "extendGap" : 1}))
##scorer2 = aminoAcidScorer("blosum62mixte", dict({"openGap" : 6, "extendGap" : 1}))
#test_SPS('BB11001', scorer1)