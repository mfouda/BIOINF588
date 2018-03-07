from blocalignments import aligne_multiple
from blocs import bloc
from seqStruct import seqStruct
from Bio import SeqIO


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

#transforme un fichier msf (alignement de référence) en dictionnaire de SeqStruct
def msftoDict (filename) :
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
                j2 += 1
                if(j2>=nb_col):
                    break

    return score/score_ref


def test_SPS(filename, scorer) :
    fasta_seqs = SeqIO.parse(open('../RV11/' + filename + '.tfa'), 'fasta')
    seqs = []
    for fs in fasta_seqs :
        seqs.append(strToSeqStruct(fs.id, str(fs.seq)))
    our_bloc = aligne_multiple(seqs, scorer)
    ref = msftoDict('../RV11/' + filename + '.msf')
    print(SPS(ref, our_bloc))

import score
scorer1 = score.aminoAcidScorer("blosum62", dict({"openGap" : 11, "extendGap" : 1}))
##scorer2 = aminoAcidScorer("blosum62mixte", dict({"openGap" : 6, "extendGap" : 1}))
test_SPS('BB11001', scorer1)

