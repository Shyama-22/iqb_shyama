#inputting penalties with score
gap_p=-1
match_sc=+2
mismatch_p=-3
#function for checking global alignments
def global_aligment(s_a,seq_b):
    rr=len(s_a)+1
    cc=len(seq_b)+1
    dp = [ [0]*(cc) for i in range(rr)]
    

    for tle in range(1,(rr)):
        dp[tle][0]=dp[tle-1][0]+gap_p
    for tle in range(1,(cc)):
        dp[0][tle]=dp[0][tle-1]+gap_p
    for i in range(1,rr):
        for j in range(1,cc):
            if(s_a[i-1]==seq_b[j-1]):
                dp[i][j]=dp[i-1][j-1]+match_sc
            else:#here we have  choices either mismatch or select a gap
                left_gap_penalty=dp[i-1][j]+gap_p
                up_gap_penalty=dp[i][j-1]+gap_p
                mismatched_penalty=dp[i-1][j-1]+mismatch_p
                dp[i][j]=max(mismatched_penalty,max(left_gap_penalty,up_gap_penalty))
    return dp          
#backtracking inorder to find optimal alignments
def optimal_s1(s_a,seq_b,dp,align_seq_a,align_seq_b,i,j,abc,maxm):
    if(i==0 and j==0):
        if(abc==maxm):
            pqr=align_seq_b[::-1]
            xyz=align_seq_a[::-1]
            print(xyz)
        
            for i in range(len(pqr)):
                if  xyz[i] != pqr[i] or xyz[i] == pqr[i] == "-":
                    print(' ', end = '')
                elif xyz[i] == pqr[i]:
                    print('|', end = '')
            print()
            print(pqr) 
            print("===========")
            return
        else:
            return         
    if (i <0 or j< 0):
        if(abc==maxm):
            print(align_seq_a[::-1]+"\n"+align_seq_b[::-1])
            print("=========")
            return

    if(i==0):
        return optimal_s1(s_a,seq_b,dp,align_seq_a+"-",align_seq_b+seq_b[j-1],i,j-1,abc+gap_p,maxm)
    if(j==0):
        return optimal_s1(s_a,seq_b,dp,align_seq_a+s_a[i-1],align_seq_b+"-",i-1,j,abc+gap_p,maxm)
    else:
        if(s_a[i-1]==seq_b[j-1]):
            optimal_s1(s_a,seq_b,dp,align_seq_a+s_a[i-1],align_seq_b+seq_b[j-1],i-1,j-1,abc+match_sc,maxm)
        else:
            optimal_s1(s_a,seq_b,dp,align_seq_a+s_a[i-1],align_seq_b+seq_b[j-1],i-1,j-1,abc+mismatch_p,maxm)

        
        optimal_s1(s_a,seq_b,dp,align_seq_a+s_a[i-1],align_seq_b+"-",i-1,j,abc+gap_p,maxm)

        optimal_s1(s_a,seq_b,dp,align_seq_a+"-",align_seq_b+seq_b[j-1],i,j-1,abc+gap_p,maxm)
#seq input
s_a="GATGCGCAG"
seq_b="GGCAGTA" 

i=len(s_a)
j=len(seq_b)
dp=global_aligment(s_a,seq_b)

print("==============================================================")
#priting all optimal alignments with score
optimal_s1(s_a,seq_b,dp,"","",i,j,0,dp[len(s_a)][len(seq_b)])
print("The alignment score for the above sequences with given penalties is:",dp[len(s_a)][len(seq_b)])
