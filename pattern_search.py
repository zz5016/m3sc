S = "ATGTTGTACCGTATCGGG"
P = "GTA"

nS = len(S)
nP = len(P)
ind = []
for i in range(nS-nP+1):
    match = True
    for j in range(nP):
        if P[j] != S[i+j]:
            match = False
            break
    if match:
        ind.append(i)
# cost = O(nS*nP)



# Problem setup

# Set pattern, P

# Read in gene sequence, S

# Map strings to base 4
def string2base4(A):
    return Abase4

def hash_eval(Xbase4):
    return Xbase10

def hash_update(HSi):
    return HSiplus1
    
# Evaluate H for pattern and first substring of S
# Pre-compute 4**nP

# Iterate through Sbase4
    # Compare HP, HSi
        # if match:
            # Compare Pbase4, Sibase4
                # if match: store index

    # Compute HSiplus1
