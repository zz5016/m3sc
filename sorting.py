import numpy as np
L = np.sort(np.random.randint(1,20,8))

n = len(L)
S = L.copy()

for i in range(1,n):
    for j in range(i-1,-1,-1):
        if S[j+1] < S[j]:
            [S[j],S[j+1]] = [S[j+1],S[j]]
        else:
            break
# cost = O(n^2)

def merge(L,R):
    nL,nR = len(L),len(R)
    indL,indR = 0,0
    M = []

    for i in range(nL+nR):
        if L[indL] < R[indR]:
            M.append(L[indL])
            indL = indL+1
            if indL == nL:
                M.extend(R[indR:])
                break
        else:
            M.append(R[indR])
            indR = indR+1
            if indR == nR:
                M.extend(L[indL:])
                break
    return M

def mergesort(X):
    n = len(X)
    if n == 1:
        return X
    else:
        L = mergesort(X[:n//2])
        R = mergesort(X[n//2:])
        return merge(L,R)
# cost = O(nlog_2(n))
