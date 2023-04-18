#function to find global alignment
def global_alignment(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    p = len(seq1)
    q = len(seq2)
    rr = len(seq1) + 1
    cc = len(seq2) + 1
    dp = []
    for i in range(rr):
        r1 = [0] * (cc)
        dp.append(r1)
    for i in range(rr):
        dp[i][0] = i * gap_penalty
    for j in range(cc):
        dp[0][j] = j * gap_penalty

    for i in range(1, rr):
        for j in range(1, cc):
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + match_score
            else:#we have 3 choices->mismatch or select a gap
                upper=dp[i-1][j]+gap_penalty
                left_align=dp[i][j-1]+gap_penalty
                diag_aln=dp[i-1][j-1]+mismatch_penalty
                counter = max(upper, left_align, diag_aln)
                dp[i][j] = counter
                
    print("Resulting DP Matrix:")
    print("===================================================")
    for row in dp:
        print(row)
    print("====================================================")

    max_score = dp[p-1][q-1]
    return max_score
sequences = {
    "seq1": "GATGCGCAG",
    "seq2": "GGCAGTA",
    "match_score": +2,
    "mismatch_penalty": -1,
    "gap_penalty": -3
}

max_score = global_alignment(
    sequences["seq1"],
    sequences["seq2"],
    sequences["match_score"],
    sequences["mismatch_penalty"],
    sequences["gap_penalty"]
)
#printing max score of alignment
print("================================================")
print("Highest score of the alignment is =", max_score)
print("====================================================")