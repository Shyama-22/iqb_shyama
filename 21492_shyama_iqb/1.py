#function for finding global alignment
def global_alignment(s_1, s_2, match_score, mismatch_p, gap_p):
    p = len(s_1)
    q = len(s_2)
    rr = len(s_1) + 1
    cc = len(s_2) + 1
    dp = []
    for i in range(rr):
        r1 = [0] * (cc)
        dp.append(r1)
    for i in range(rr):
        dp[i][0] = i * gap_p
    for j in range(cc):
        dp[0][j] = j * gap_p
#here we are computing the scores of upper,left and diagonal alignments
    def upper(i, j):
        return dp[i-1][j]+gap_p
    
    def left_align(i, j):
        return dp[i][j-1]+gap_p
    
    def diag_aln(i, j):
        if s_1[i - 1] == s_2[j - 1]:
            return dp[i - 1][j - 1] + match_score
        else:
            return dp[i-1][j-1]+mismatch_p

    for i in range(1, rr):
        for j in range(1, cc):
            dp[i][j] = max(upper(i, j), left_align(i, j), diag_aln(i, j))

    #printing matrix
    print("Resulting DP Matrix:")
    print("===================================================")
    for i in range(rr):
        for j in range(cc):
            print(dp[i][j], end=" ")
        print()
    print("=====================================================")

    max_score = dp[p][q]
    return max_score

#inputting sequence with resp scores
sequences = {
    "s_1": "GATGCGCAG",
    "s_2": "GGCAGTA",
    "match_score": +2,
    "mismatch_p": -3,
    "gap_p": -1
}

max_score = global_alignment(
    sequences["s_1"],
    sequences["s_2"],
    sequences["match_score"],
    sequences["mismatch_p"],
    sequences["gap_p"]
)
#printing the highest score of alignment
print("Highest score of the alignment is =", max_score)


