"""M345SC Homework 1, part 2
Name: Ziyao Zhang
CID: 01181949
"""

def nsearch(L,P,target):
    """
    Input:
    L: list containing *N* sub-lists of length M. Each sub-list contains M
       numbers (floats or ints), and the first P elements of each sub-list can
       be assumed to have been sorted in ascending order (assume that P<M).
       L[i][:p] contains the sorted elements in the i+1th sub-list of L.
    P: the first P elements in each sub-list of L are assumed to be sorted in
       ascending order target: the number to be searched for in L.

    Output:
    Lout: a list consisting of Q 2-element sub-lists where Q is the number of
          times target occurs in L. Each sub-list should contain 1) the index of
          the sublist of L where target was found and 2) the index within the
          sublist where the target was found. So, Lout = [[0,5],[0,6],[1,3]]
          indicates that the target can be found at L[0][5],L[0][6],L[1][3]. If
          target is not found in L, simply return an empty list (as in the code
          below).

    Discussion:
    - Algorithm
        The function searches for targets in a combination of binary and linear
    fashion. For each sequence in the list, it searches binarily for the first
    and the last target in the sorted part one after another, and searches
    linearly for the rest targets in the unsorted part of the sequence. The
    binary search is terminated when some target found is acsertained to be
    adjacent to a non-target element. At last, append all locations in between
    the first and the last target found by binary search to Lout, alongside the
    ones in the unsorted part found by linear search.

    - Running time analysis
        The computational cost is O(log_2(P)) for each binary search and O(M-P)
    for each linear search. To record all locations of the targets, the function
    has to run through all found locations. This is potentially a O(M) running
    time in the worst case scenario (for example a list consists of repeated
    targets only).

    In conclusion, the algorithm above generally has an O(MN) time complexity as
    the process is repeated for every sublist. However the cost can be vastly
    reduced in some 'well-structured' cases. It happens when there are few
    targets showing up in the list. This apparently saves the cost of recording
    the locations, which according to previous analysis, is the major cause of
    the leading order term. It also happens for large P/M, i.e. a large
    proportion of the sequence is sorted, in which case the cost for binary
    search dominates in running time.

    This algorithm can be considered efficient as it makes use of binary search
    for the sorted part of the list. Compared to looking for a target and then
    linearly searching in its neighbourhood, it outperforms linear search when
    it comes to the worst case (O(log_2(n)) vs O(n)) and saves largely the cost
    for the cases which have large proportion of reduplication of targets.
    """

    # setup
    Lout = []
    for i in range(len(L)):
        subL = L[i]
        istart = 0
        iend = P-1
        left = -100

        # binary search for the first target in the sorted sub-sequence
        while istart <= iend:
            imid = (istart+iend)//2
            if subL[imid] < target:
                if subL[imid+1] == target:
                    left = imid+1
                    break
                if subL[imid+1] < target:
                    istart = imid+1
                else:
                    break
            else:
                if imid == 0:
                    if subL[imid] == target:
                        left = imid
                        break
                else:
                    iend = imid-1

        # binary search for the last target in the sorted sub-sequence
        if left != -100:
            istart = left+1
            iend = P-1
            right = left
            while istart <= iend:
                imid = (istart+iend)//2
                if subL[imid] > target:
                    if subL[imid-1] == target:
                        right = imid
                        break
                    else:
                        iend = imid-1
                else:
                    if imid == P-1:
                        if subL[imid] == target:
                            right = imid+1
                            break
                    else:
                        istart = imid+1

            # store locations
            for n in range(left,right):
                Lout.append([i,n])

        # linear search for targets in the unsorted sub-sequence
        for ind,x in enumerate(subL[P:]):
            if x == target:
                Lout.append([i,ind+P])

    return Lout


if __name__ == '__main__':
    L = [[2, 2, 2, 2, 1, 4, 6, 2, 4, 2, 1], [1, 2, 2, 4, 1, 4, 6, 2, 4, 3, 1]]
    print(nsearch(L,4,2))
    print(L)
