import numpy as np
L = np.sort(np.random.randint(1,20,8))
print("L =",L)

x = 7
ind = "not found"
for i,l in enumerate(L):
    if x == l:
        ind = i
        break
# cost = O(n)

def search(L,x):
    """
    Use binary search to find targe x in sorted array L.
    """

    assert type(x) is int or type(x) is float, "error, x must be numeric"
    # check if L is sorted...

    ind = "not found"
    n = len(L)
    istart = 0
    iend = n-1
    while istart <= iend:
        imid = (iend+istart)//2
        if x == L(imid):
            ind = imid
            break
        elif x < L(imid):
            iend = imid-1
        else:
            istart = imid+1
    return ind
# cost = O(log_2(n))
