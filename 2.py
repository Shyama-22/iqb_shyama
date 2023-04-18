keys_map =['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
a1_v =[1.45,0.77,0.98,1.53,1.12,0.53,1.24,1,1.07,1.34,1.2,0.73,0.59,1.17,0.79,0.79,0.82,1.14,1.14,0.61]
b1_v=[0.97,1.30,0.80,0.26,1.28,0.81,0.71,1.6,0.74,1.22,1.67,0.65,0.62,1.23,0.9,0.72,1.2,1.65,1.19,1.29]
alpha = {keys_map[i]: a1_v[i] for i in range(len(keys_map))}
beta = {keys_map[i]: b1_v[i] for i in range(len(keys_map))}
#a1_sc(seq): This function takes a protein sequence as input and returns the sum of the alpha values 
# of each amino acid in the sequence. The alpha values are stored in a dictionary called alpha 
# where each amino acid is associated with its corresponding alpha value.
def a1_sc(seq):  
    return sum(alpha[i] for i in seq)
#b1_sc(seq): This function takes a protein sequence as input and returns the sum of the beta 
# values of each amino acid in the sequence. The beta values are stored in a dictionary called beta
# where each amino acid is associated with its corresponding beta value.
def b1_sc(seq):  
    return sum(beta[i] for i in seq)
#The function "h_p_a_v(i, p_val)" takes two arguments i and p_val and returns the number of consecutive 
# amino acids from position i to i+p_val that have an alpha score of 1 or more.
def h_p_a_v(i, p_val):
    return sum(1 for j in range(i, i + p_val) if a1_sc(seq[j]) >= 1)
#he function "h_p_b_v(i, p_val)" takes two arguments i and p_val and returns the number of consecutive 
# amino acids from position i to i+p_val that have a beta score of 1 or more.
def h_p_b_v(i, p_val):
    return sum(1 for j in range(i, i + p_val) if b1_sc(seq[j]) >= 1)  
#The function "a_extending(arr, p1, p2, length_seq)" extends a list arr of H (helix) characters on both ends 
# of the protein sequence based on the alpha scores of the amino acids in the sequence.
def a_extending(array,p1,p2,length_seq):
    while p1 < length_seq and a1_sc(seq[p1-3:p1+1]) >= 4:
        array[p1] = 'H'
        p1 += 1
    while p2 >= 0 and a1_sc(seq[p2:p2+4]) >= 4:
        array[p2] = 'H'
        p2 -= 1
# function b_extending(arr, p1, p2, length_seq)" extends a list arr of S (sheet) characters on both ends
# of the protein sequence based on the beta scores of the amino acids in the sequence.
def b_extending(array,p1,p2,length_seq):
    while p1 < length_seq and b1_sc(seq[p1-3:p1+1]) >= 4:
        array[p1] = 'S'
        p1 += 1
    while p2 >= 0 and b1_sc(seq[p2:p2+4]) >= 4:
        array[p2] = 'S'
        p2 -= 1

#The function "alfa__helix(seq)" predicts helix regions in a protein sequence by searching for 6
# consecutive amino acids with alpha score >= 1, and then extends these regions using the "extend_a" function.
def alfa__helix(seq):
    if_helixing = ['_' for _ in range(len(seq))]
    
    if len(seq) <= 5:
        raise ValueError("Sequence length should be greater than 5.")
    
    for i in range(len(seq)-5):
        n = sum([1 for j in range(i, i+6) if a1_sc(seq[j]) >= 1])
        if n >= 4:
            for j in range(i, i+6):
                if if_helixing[j] != 'H':
                    if_helixing[j] = 'H'
            p1, p2 = i+6, i-1
            a_extending(if_helixing, p1, p2, len(seq))

    return "".join(if_helixing)
#The function "beta_sheet(seq)" predicts sheet regions in a protein sequence by searching for 5 
# consecutive amino acids with beta score >= 1, and then extends these regions using the "extend_b" function.
def beta_sheet(seq): 
    sheet = ["_" for _ in range(len(seq))]
    for i in range(len(seq) - 4):
        if h_p_b_v(i, 5) >= 3:
            for j in range(i, i+5):
                if sheet[j] != 'S':
                    sheet[j] = 'S'
            b_extending(sheet, i+5, i-1, len(seq))
    return "".join(sheet)
#The function "funcs(seq1, seq2, seq, res)" compares the predicted helix and sheet regions of two aligned protein sequences
# seq1 and seq2, and returns the resulting alignment in res. 
def funcs(sequence1,sequence2,seq,resulting):
    tle=len(seq)
    i=0
    while i < tle:
        if (sequence1[i] == 'H' and sequence2[i] == '_') or (sequence1[i] == '_' and sequence2[i] == 'H'):
            resulting[i] = 'H'
            i += 1
        elif (sequence1[i] == '_' and sequence2[i] == 'S') or (sequence1[i] == 'S' and sequence2[i] == '_') :
            resulting[i] = 'S'
            i += 1
        elif sequence1[i] == '_' and sequence2[i] == '_':
            resulting[i] = '*'
            i += 1
        else:
            n = 0
            while (sequence1[i] == 'H' and sequence2[i] == 'S') or (sequence1[i] == 'S' and sequence2[i] == 'H'):
                n += 1
                if i < tle-1:
                    i += 1
                else:
                    break
            if i == tle-1:
                i += 1 
            p1 = a1_sc(seq[i-n:i])
            p2 = b1_sc(seq[i-n:i])
            if p1 > p2:
                for k in range(i-n,i):
                    resulting[k] = 'H'
            else:
                for k in range(i-n,i):
                    resulting[k] = 'S'
    return resulting              
    
def conflict_case(seq1, seq2, seq):
    result = ["" for i in range(len(seq))]
    funcs(seq1, seq2, seq, result)
    return "".join(result)
#here we r inputting our sequence 
seq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF" 
s1 = alfa__helix(seq)  
s2 = beta_sheet(seq)  
resulting_sq = conflict_case(s1, s2, seq)  
print("----------------------------------------------")
print("our sequence is="+seq)
#finalyy printing our resulting seq
print("our resulting seq is="+resulting_sq)
