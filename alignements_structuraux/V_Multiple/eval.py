from blocalignments import aligne_multiple
from blocs import bloc



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


def parse_msf (filename) :
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
    while (i + nb_seqs - 1) < n :
        for k in range (nb_seqs) :
            line = lines[i+k]
            l = len(line)
            j = 0
            while j<l and not line[j]==' ' :
                j += 1
            while j<l and line[j]==' ' :
                j += 1
            seqs[k] += line[j:l-1]
        i += nb_seqs + 2
    
    for k in range (nb_seqs) :
        gap = False
        seqs[k] = seqs[k].replace(' ','')
        seqs[k] = seqs[k].replace('.','-')
        for c in range(len(seqs[k])) :
            if seqs[k][c] == '-' :
                if not gap :
                    seqs[k] = seqs[k][:c] + '+' + seqs[k][c+1:]
                    gap = True
            else :
                gap = False
    return seqs


def score_C2 (filename, d, e) :
    seqs_input = parse_tfa((filename + '.tfa'))
    bloc_result = aligne_multiple(seqs_input, d, e)
    aligne_ref = parse_msf(filename + '.msf')
    
    n = bloc_result.getNbSeqs()
    m = len(bloc_result.getSeqs()[0])
    score_C = 0
    score_C2 = 0
    
#    fichier = open(filename + '_score_C2.txt', 'w')
#    carac_per_line = 70
#    lines = ['']*(2 + m//carac_per_line)*(2*n + 8)
#    line = 0
#    carac = 0
    
    for j in range (m) :
        col_nous, col_ref = ordonne_col(j, n, aligne_ref, bloc_result)
        score_col = 0
        
#        if carac == carac_per_line :
#            line += 2*n + 6
#            carac = 0
        
        for i in range (n) :
#            lines[line + i] += col_nous[i] + ' '
#            lines[line + i + n + 1] += col_ref[i] + ' '
            score_col += (col_nous[i] == col_ref[i])
        
#        lines[line + 2*n + 2] += str(score_col) + ' '
#        carac += 1           
        score_C2 += score_col
        score_C += (score_col == n)
    
#    lines[line + 2*n + 5] += 'score_C2 = ' + str(score_C2) + ' / ' + str(m*n) + ' = ' + str(score_C2/(m*n))
#    lines[line + 2*n + 4] += 'score_C = ' + str(score_C) + ' / ' + str(m) + ' = ' + str(score_C/m)
#    
#    for x in lines :
#        fichier.write(x + '\n')
#    fichier.close()
#    
    return score_C/m, score_C2/(m*n)


m = 0
for i  in range(38) :
    zero = ''
    if i<9 :
        zero = '0'
    if i != 4 :
        s = score_C2('../RV11/BBS110' + zero + str(i+1), 6, 1)[0]
        print(s)
        m += s
print(m/37)