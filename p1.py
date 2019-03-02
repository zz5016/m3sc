"""M345SC Homework 1, part 1
Name: Ziyao Zhang
CID: 01181949
"""

from time import time

def ksearch(S,k,f,x):
    """
    Search for frequently-occurring k-mers within a DNA sequence and find the
    number of point-x mutations for each frequently occurring k-mer.

    Input:
    S: a string consisting of A,C,G, and Ts.
    k: the size of k-mer to search for (an integer).
    f: frequency parameter -- the search should identify k-mers which occur at
       least f times in S.
    x: the location within each frequently-occurring k-mer where point-x
       mutations will differ.

    Output:
    L1: list containing the strings corresponding to the frequently-occurring
        k-mers.
    L2: list containing the locations of the frequent k-mers.
    L3: list containing the number of point-x mutations for each k-mer in L1.

    Discussion:
    - Algorithm
        The function runs through a given DNA sequence and creates a dictionary
    that records patterns of required length as keys and stores their locations
    in lists as values.(D) The function extracts so-called frequently-occurring
    k-mers from the dictionary by checking if their location lists are no less
    than the desired frequency.(L1) Obtain simliarly a list of location lists
    for frequently-occurring patterns.(L2)

    Furthermore, based on the existent dictionary, the function constructs a new
    dictionary where keys are set to be sub-patterns with the nucleotide at
    point-x removed, and values are set to be the numbers of occurrence of those
    sub-patterns. They are computed by summing up the lengths of location lists
    of their correponding original patterns in the previous dictionary. This
    figure also represents the total frequency of a pattern and all its point-x
    mutations.(d)

    Repeat the process of generating sub-patterns for frequently-occurring
    k-mers and read their values from the new dictionary. The function works out
    the number of mutations in the sequence of a frequently-occurring pattern by
    subtracting its own frequency from the read value.(L3)

    - Running time analysis
        The first part of above algorithm runs through almost every single
    position in the sequence (assume k-mers are far shorter than the sequence).
    At each position it draws k elements from the sequence, makes a dictionary
    search (O(1)) and either creates a new key or appends current position to
    the location lists. Therefore, the construction of the first dictionary
    requires O(len(S)*k) running time (up to some constant). There is no much
    difference between best and worst cases.

    The construcion of L1 requires the function to run through all keys in the
    dictionary. With each key it calculates a length and makes a comparison.
    Therefore its computational cost is strictly less than some multiple of
    len(S). This can be largely reduced if the sequence contains a vast amount
    of reduplicative k-mers. (This happens for small k, in which case the time
    complexity drops to O(4^k), as 4^k is the possible number of distinct
    k-mers.)

    L2 basically contains a copy of positions (potentially all positions in the
    sequence) from the dictionary and all elements in L1 are called as positions
    are added into L2. As a result, it brings another O(len(S)) running time.
    However, as it is impossible to have all k-mers frequently-occurring, the
    cost can be massively saved, especially for large k or large f.

    To construct the second dictionary, the function goes for every key in the
    first dictionary, and then extracts k-1 elements, makes a dictionary search,
    calculates a length and updates a value. The function further calls again
    all elements in L1 and similarly each call is followed by O(k) operations.
    In conclusion, this part has an O(len(S)*k) (or O(k*4^k) for small k) time
    complexity in total. As stated above, it would go far less that an actual
    O(len(S)*k) running time if the sequence contains a large proportion of
    repeated k-mers.

    - Summary
        To sum up, the function runs through the whole sequence once and both
    first dictionary and L1 twice. The leading-order running time is obviously
    len(S)*k. It is as anticipated that this term plays the dominant part,
    because it is caused by running the whole sequence and is a necessary cost
    that has almost little chance to reduce.

    This algorithm can be considered as efficient since it takes advantage of
    dictionaries. Although the drawback is that it may require numerous memory
    space in some extreme cases (my own laptop does not have a large RAM and it
    runs into a memory error when it comes to a very long sequence), it largely
    shortens the searching time cost.

    Moreover, there are many different scenarios can that possibly save some
    terms (other than the leading-term though, unfortunately) of running time,
    both for large k or small k and large f, which suggests that the algorithm
    works well in general. (Experiments showed that it generates the result
    correctly and efficiently, and performs extremely fast for recurrent DNA
    sequences with small k.)
    """

    # dictionary I
    D = {}
    nS = len(S)
    for i in range(nS-k+1):
        P = S[i:(i+k)]
        if P in D:
            D[P].append(i)
        else:
            D[P] = [i]

    # L1
    L1 = [P for P in D if len(D[P]) > f-1]

    # L2
    L2 = [D[P] for P in L1]

    # dictionary II
    d = {}
    for P in D:
        M = P[:x] + P[x+1:]
        if M in d:
            d[M] += len(D[P])
        else:
            d[M] = len(D[P])

    # L3
    L3 = []
    for P in L1:
        M = P[:x]+P[x+1:]
        L3.append(d[M]-len(D[P]))

    return L1,L2,L3


if __name__ == '__main__':
    S = 'CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    k = 3
    x = 2
    f = 2
    print(ksearch(S,k,f,x))

    S = 'CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'*10**5
    k = 4
    x = 0
    f = 2
    t1 = time()
    ksearch(S,k,f,x)
    t2 = time()
    print(t2-t1)
