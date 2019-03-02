"""M345SC Homework 2, part 1
Name: Ziyao Zhang
CID: 01181949
"""

def scheduler(L):
    """
    Question 1.1
    Schedule tasks using dependency list provided as input.

    Input:
    L: Dependency list for tasks. L contains N sub-lists, and L[i] is a sub-
       list containing integers (the sub-list my also be empty). An integer, j,
       in this sub-list indicates that task j must be completed before task i
       can be started.

    Output:
    S: A list of integers corresponding to the schedule of tasks. L[i]
       indicates the day on which task i should be carried out. Days are
       numbered starting from 0.

    Discussion:
    - Algorithm:
        It can be proved that the earliest day on which a task can be possibly
    carried out is the maximal distance between the task and all the independent
    tasks. This can be computed recursively since the distance of a task is
    simply the maximum of the distances of its pre-tasks plus one. The inner
    function d computes the desired distance in such a way and meanwhile stores
    every distance called due to recursion (in the list S). When calling any
    explored distance, the function reads the distance directly from S instead
    of using recursion again, so that repetitions can be avoided. Lastly let the
    function d run through all tasks.

    - Running Time and Efficiency:
        In the following context, N = len(L) refers to the number of tasks and
    M = sum of len(sublists of L) refers to the number of dependence between
    tasks.

    1) The function runs through all tasks and for each task it either computes
       the distances for its unexplored pre-tasks, or reads the distances from S
       for the pre-tasks that are already explored. This process allows all the
       links between the tasks to be covered for constant times and hence in
       total has a computational cost of O(M). For the similar reason, looking
       for the maximum among the distances inherited from pre-tasks also
       requires O(M) running time.
    2) Both assigning distances in their correponding position in S and calling
       the inner functions in the loop bring another O(N) time complexity as
       there are N tasks to be completed.
    3) Overall this algorithm has a O(M+N) time complexity in general. The best
       and the worst case do not make a big difference.

    The algorithm can be considered as efficient since its running time depends
    linearly on N and M, which is a necessary cost and has no reason to reduce.
    """
    n = len(L)
    S = [-1]*n

    # maximal distance of a task from initiators
    def d(i):
        if S[i] != -1:
            pass
        elif len(L[i]) == 0:
            S[i] = 0
        else:
            # distances inherited from pre-tasks
            dist = []
            dist.extend([d(j)+1 for j in L[i] if S[j]==-1])
            dist.extend([S[j]+1 for j in L[i] if S[j]!=-1])
            # store distance in S
            S[i] = max(dist)
        return S[i]

    # compute distance for all tasks
    for i in range(n):
        d(i)

    return S



def findPath(A,a0,amin,J1,J2):
    """
    Question 1.2.i
    Search for feasible path for successful propagation of signal from node J1
    to J2.

    Input:
    A: Adjacency list for graph. A[i] is a sub-list containing two-element
       tuples (the sub-list my also be empty) of the form (j,Lij). The integer,
       j, indicates that there is a link between nodes i and j and Lij is the
       loss parameter for the link.
    a0: Initial amplitude of signal at node J1
    amin: If a>=amin when the signal reaches a junction, it is boosted to a0.
          Otherwise, the signal is discarded and has not successfully reached
          the junction.
    J1: Signal starts at node J1 with amplitude, a0.
    J2: Function should determine if the signal can successfully reach node J2
        from node J1.

    Output:
    L: A list of integers corresponding to a feasible path from J1 to J2.

    Discussion:
    - Algorithm:
        Two junctions are considered as connected if a signal of amplitude a0
    that is weakened via direct propagation between them still exceeds the
    threshold. The function searches aggresively and recursively for the target
    J2 from the source J1 in a similar way as depth first search.

    The funciton updates the status of unexplored neighbours as explored while
    searching. It add the neighbours one by one to the end of the path, or, if
    the path meets a dead end (no unexplored junctions are connected to the
    current location), it removes the current location from the path. Once a
    junction is added, the above process is repeated immediately on the most
    recently added location via recursion. The function stops further recursions
    and breaks from all unfinished loops and returns the path when J2 is found,
    or returns an empty list if no more junctions can be explored.

    - Running Time and Efficiency:
        In the following context, N = len(A) refers to the number of junctions
    and M = sum of len(sublists of A)/2 refers to the number of connections
    between junctions.

    Each time the inner function findJ2 is called:
    1) The function checks if a path to J2 is already found by the indicator I
       for each iteration. It runs through all neighbours of current junction by
       the for loop and only the unexplored neighbours are added to the path and
       marked as explored at different stages of the process (as each iteration
       in the for loop triggers the recursion to happen). This in total requires
       a running time proportional to the number of connections.
    2) It removes the last element in the path unless J2 is found, which costs
       constant time (O(1)).
    3) The use of recursion causes the above computations to be carried out on
       possibly all junctions that are in the same component as J1 (and once each
       because it runs only if a junction is unexplored). In other word, the
       inner function needs to call itself at most M times to cover all such
       junctions. As a result, the cost in total depends linearly on M both for
       (1) and (2).

    In conclusion, the overall computational cost is O(M). Compared to the
    original method of depth-first search, the above algorithm avoids storing
    and passing the path for every node, which requires both large memory and
    slows down the computation as producing a copy of path brings O(M^2)
    operations. Therefore in short, this algorithm is efficient in general
    (though in some cases where J2 is very close to J1, breadth-first search
    would performs better).
    """
    Lmin = amin/a0
    L = [J1]

    # special case
    if J1 == J2:
        return L

    E = {J1:1}  # explored junctions
    I = {J2:False}  # indicator

    # aggressive search for J2 from i
    def findJ2(i):
        for j,l in A[i]:
            # break and avoid further recursions if L found
            if I[J2]:
                break
            # search through neighbours
            if j not in E and l >= Lmin:
                L.append(j)
                E[j] = 1
                # J2 found
                if j == J2:
                    I[J2] = True
                    break
                # find path to J2 by recursion
                findJ2(j)
        # backtrack if current path meets a dead end
        if not I[J2]:
            L.pop()

    findJ2(J1)
    return L



def a0min(A,amin,J1,J2):
    """
    Question 1.2.ii
    Find minimum initial amplitude needed for signal to be able to successfully
    propagate from node J1 to J2 in network (defined by adjacency list, A).

    Input:
    A: Adjacency list for graph. A[i] is a sub-list containing two-element
       tuples (the sub-list my also be empty) of the form (j,Lij). The integer,
       j, indicates that there is a link between nodes i and j and Lij is the
       loss parameter for the link.
    amin: Threshold for signal boost. If a>=amin when the signal reaches a
          junction, it is boosted to a0. Otherwise, the signal is discarded and
          has not successfully reached the junction.
    J1: Signal starts at node J1 with amplitude, a0.
    J2: Function should determine min(a0) needed so the signal can successfully
        reach node J2 from node J1.

    Output:
    (a0min,L) a two element tuple containing:
    a0min: minimum initial amplitude needed for signal to successfully reach J2
           from J1.
    L: A list of integers corresponding to a feasible path from J1 to J2 with
       a0=a0min. If no feasible path exists for any a0, return output as shown
       below.

    Discussion:
    - Algorithm:
        The problem basically looks for the maximal Lmin such that there exists
    a path from J1 to J2 in which the minimum of Lij's is no less than Lmin. The
    idea of the following algorithm refers to Dijkstra's algorithm, except for
    the prioity - it chooses to move to an unexplored neighbour along the link
    with maximal weight (loss parameter). This algoritm guarantees to find the
    maximal 'minimal weight' among all paths from J1 to J2, which is equal to
    the desired Lmin.

    The function stores the path with maximal 'minimal weight' from J1 to every
    newly found junction in a dictionary and their minimal weights in another
    dictionary. The path and the weight get updated if an alternative path to
    the junction with larger 'minimal weight' is found. A junction is marked
    explored only if the path to the junction has the maximal 'minimal weight'
    among all unexplored junctions in the dictionary. It is then removed from
    both dictionaries. The function repeats above process and terminates when J2
    is marked as explored.

    - Running Time and Efficiency:
        In the following context, N = len(A) refers to the number of junctions
    and M = sum of len(sublists of A)/2 refers to the number of connections
    between junctions.

    For each iteration in the while loop, non-constant costs are as follow:
    1) It takes the key of maximal value from the weight dictionary, which
       involves O(N) operations as the dictionary contains N items at most.
    2) The function updates paths by duplicating the path to the current
       junction i and passing it its unexplored neighbours. This causes O(Mj*N)
       time complexity, where Mi is the number of connections at junction j.

    The dictionary contains at most N items and for each iteration an item gets
    removed. As a result the while loop runs N times in the worst case.
    Therefore the algorithm requires O(N*N+sum(Mi)*N) = O(N^2+M*N) time cost.
    (In the case where M<<N, the time complexity is reduced to O(M^2) as in this
    case both the dictionary and paths contain at most M elements instead.
    Therefore it makes more sense to replace N by M.)
    """
    # special case
    if J1 == J2:
        a0min,L = amin,[J1]
        output = a0min,L
        return output

    # function that obtains key with largest value from dictionary
    def maxkey(d):
        vmax = -1
        for k,v in d.items():
            if v > vmax:
                vmax = v
                kmax = k
        return kmax

    a0min,L = -1,[]
    E = {}  # explored junctions
    path = {}  # paths
    w = {}  # weights

    path[J1] = [J1]
    w[J1] = 1

    while len(w) != 0:
        i = maxkey(w)
        E[i] = 1
        # arrive J2 via a path with maximal weight
        if i == J2:
            a0min,L = amin/w[i],path[i]
            break
        # search through and update information for neighbours of current location
        for j,l in A[i]:
            if j not in E and (j not in w or w[j] < min(l,w[i])):
                # add or update path and weight
                w[j] = min(l,w[i])
                path[j] = path[i] + [j]
        # delete explored junction
        del w[i]
        del path[i]

    output = a0min,L
    return output



if __name__ == '__main__':
    # test
    L = [[],[0],[0],[1],[0,1,2,5],[2],[0,1,4],[6],[3,7]]
    print(scheduler(L))

    A = [[(1,0.7),(2,0.4),(3,0.7)],
         [(0,0.7)],
         [(0,0.4),(3,0.6)],
         [(0,0.7),(2,0.6)]]
    B = [[(1,0.8),(2,0.8),(3,0.9)],
         [(0,0.8),(4,0.5)],
         [(0,0.8),(4,0.7)],
         [(0,0.9),(4,0.6)],
         [(1,0.5),(2,0.7),(3,0.6),(5,0.7),(6,0.6),(9,0.3)],
         [(4,0.7),(6,0.7)],
         [(4,0.6),(5,0.7),(7,0.4),(8,0.5)],
         [(6,0.4),(10,0.2)],
         [(6,0.5),(10,0.1)],
         [(4,0.3),(10,0.9)],
         [(7,0.2),(8,0.1),(9,0.9)],
         []]
    print(findPath(A,10,6,0,2))
    print(findPath(B,10,3,0,10))
    print(findPath(B,10,1,0,10))
    print(findPath(B,10,8,0,11))

    print(a0min(B,0.3,0,10))
    print(a0min(B,0.7,0,4))
    print(a0min(B,7,0,6))
    print(a0min(B,0.3,0,11))
