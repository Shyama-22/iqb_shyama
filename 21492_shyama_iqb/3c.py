#here we are inputting the penalties
gap_p=-3
match_s=+2
mismatch_p=-1
#max in 2d array
def getmaxindp(dp):
    max=0
    for i in range(0,len(dp)):
        for j in range(0,len(dp[0])):
            if dp[i][j]>=max:
                max=dp[i][j]
    return max
#finding max index
def getmaxdp(dp,max):    
    ab=0
    cd=0
    lcs=[] 
    for i in range(0,len(dp)):
        for j in range(0,len(dp[0])):
            if dp[i][j]==max:
                max=dp[i][j]
                ab=i
                cd=j
                lcs.append([ab,cd])    
    return lcs
#function for checking global alignment
def local_aligment(seq_a,seq_b):
    rr=len(seq_a)+1
    cc=len(seq_b)+1
  
    dp = [ [0]*(cc) for i in range(rr)]
    
    for i in range(1,rr):
        for j in range(1,cc):
            if(seq_a[i-1]==seq_b[j-1]):
                dp[i][j]=max(dp[i-1][j-1]+match_s,0)
            else:# we have 3 choices
                left_gap_penalty=max(dp[i-1][j]+gap_p,0)
                up_gap_penalty=max(dp[i][j-1]+gap_p,0)
                mismatched_penalty=max(dp[i-1][j-1]+mismatch_p,0)
                dp[i][j]=max(mismatched_penalty,max(left_gap_penalty,up_gap_penalty))               
    return dp  
#backtracking
def optimal_s(seq_a,seq_b,dp,align_seq_a,align_seq_b,i,j,xyz,maxm):
    
    if(dp[i][j]==0):
        if(xyz==maxm):
            srm=align_seq_b[::-1]
            uno=align_seq_a[::-1]
            print(uno)
            for i in range(len(srm)):
                if  uno[i] != srm[i] or uno[i] == srm[i] == "-":
                    print(' ', end = '')
                elif uno[i] == srm[i]:
                    print('|', end = '')
            print()
            print(srm) 
            print("===========")
            return
    if(i==0 and j==0):
        if(xyz==maxm):
            srm=align_seq_b[::-1]
            uno=align_seq_a[::-1]
            print(uno)
            for i in range(len(srm)):
                if  uno[i] != srm[i] or uno[i] == srm[i] == "-":
                    print(' ', end = '')
                elif uno[i] == srm[i]:
                    print('|', end = '')
            print()
            print(srm) 
            print("===========")
            return        
    if (i <0 or j< 0):
        if(xyz==maxm):
            print(align_seq_a[::-1]+"\n"+align_seq_b[::-1])
            print("=========")
            return
        else:
            return    

    if(i==0):
        return optimal_s(seq_a,seq_b,dp,align_seq_a+"-",align_seq_b+seq_b[j-1],i,j-1,xyz+gap_p,maxm)
    if(j==0):
        return optimal_s(seq_a,seq_b,dp,align_seq_a+seq_a[i-1],align_seq_b+"-",i-1,j,xyz+gap_p,maxm)
    else:
        if(seq_a[i-1]==seq_b[j-1]):
            optimal_s(seq_a,seq_b,dp,align_seq_a+seq_a[i-1],align_seq_b+seq_b[j-1],i-1,j-1,xyz+match_s,maxm)
        else:
            optimal_s(seq_a,seq_b,dp,align_seq_a+seq_a[i-1],align_seq_b+seq_b[j-1],i-1,j-1,xyz+mismatch_p,maxm)
        optimal_s(seq_a,seq_b,dp,align_seq_a+seq_a[i-1],align_seq_b+"-",i-1,j,xyz+gap_p,maxm)

        optimal_s(seq_a,seq_b,dp,align_seq_a+"-",align_seq_b+seq_b[j-1],i,j-1,xyz+gap_p,maxm)
#inputting sequence
S1 = "GATGCGCAG"
S2 = "GGCAGTA"
dp=local_aligment(S1,S2)
#used for finding max
tle=getmaxdp(dp,getmaxindp(dp))
for ss in tle:
    maxm=dp[ss[0]][ss[1]]
    print("The score for  sequences  is:",maxm) 
    optimal_s(S1,S2,dp,"","",ss[0],ss[1],0,maxm)
