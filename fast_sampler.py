# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak, and Strahil Ristov

import time
import sys
from collections import defaultdict
from tkinter import messagebox


def calc_dists(x, E, dist, parent=-1, cur_dist=0):
    """
    Takes a tree of size N, node x, and a vector dist of size N.
    Fills dist with numbers such that dist[i] contains distance between
    nodes x and i. 
    Time complexity: O(N).
    """
    dist[x] = cur_dist
    for y in E[x]:
        if y != parent:
            calc_dists(y, E, dist, x, cur_dist + 1)


def most_distant_pair(E, C):
    """
    Takes a tree of size N and returns two most distant nodes.
    Only considers pairs where both members are in C.
    Time complexity: O(N).
    """
    N = len(E)
    dist = N * [0]
    calc_dists(list(C)[0], E, dist)
    first = max(C, key=lambda x: dist[x])
    calc_dists(first, E, dist)
    second = max(C, key=lambda x: dist[x])
    return (first, second)


def greedy_solve(E, C, K):
    """
    Takes a tree of size N, a set of nodes C and an integer K. 
    Using suboptimal greedy algorithm finds a subset of C of size K, 
    such that sum of their pair-wise distances in the tree is maximized.
    Returns a pair (vector of chosen K nodes, sum of pair-wise distances).
    Assumes 2 <= K <= N.
    Time complexity: O(NK).
    Memory complexity: O(N)
    """
    N = len(E)
    K = min(K, len(C))

    # is_chosen[i] is true iff node i is currently chosen
    is_chosen = N * [False]
    # dists[i] is sum of distances from node i to currently chosen nodes
    dists = N * [0]
    # total_dists is sum of pair-wise distances of chosen nodes
    total_dists = 0

    dist_from_x = N * [0]

    def choose(x, is_chosen, dists, dist_from_x, E):  # make node x chosen
        assert (not is_chosen[x])
        ret = dists[x]
        is_chosen[x] = True
        calc_dists(x, E, dist_from_x)
        for i in range(N):
            dists[i] += dist_from_x[i]
        return ret

    first, second = most_distant_pair(E, C)
    total_dists += choose(first, is_chosen, dists, dist_from_x, E)
    total_dists += choose(second, is_chosen, dists, dist_from_x, E)

    counter = 2
    while counter < K:
        best = None
        for i in C:
            if (not is_chosen[i]) and (best is None or dists[i] > dists[best]):
                best = i
        total_dists += choose(best, is_chosen, dists, dist_from_x, E)
        counter += 1

    chosen = filter(lambda x: is_chosen[x], range(N))
    return (chosen, total_dists)


def dp_solve(E, C, K):
    """
    Takes a tree of size N, a set of nodes C and an integer K. 
    Using optimal DP algorithm finds a subset of C of size K, 
    such that sum of their pair-wise distances in the tree is maximized.
    Returns a pair (vector of chosen K nodes, sum of pair-wise distances).
    Assumes 2 <= K <= N.
    Time complexity: O(NK^2).
    Memory complexity: O(NK)
    """
    N = len(E)
    K = min(K, len(C))
    inf = int(1e100)

    # dp and reconstruction matrix
    f = N * [0]
    r = N * [0]

    # helper vectors
    nf = K * [0]

    def calc_dp(x, parent):
        M = 0  # number of children
        for y in E[x]:
            if y != parent:
                calc_dp(y, x)
                M += 1

        # (M+1)x(K+1) matrix initialized with -oo
        f[x] = [(K+1)*[-inf] for i in range(M + 1)]
        # (M+1)x(K+1) matrix initalized with -1
        r[x] = [(K+1)*[-1] for i in range(M + 1)]

        f[x][0][0] = 0
        f[x][0][1] = 0 if x in C else -inf
#        f[x][0][1] = 0
        idx = 0
        for y in E[x]:
            if y != parent:
                for i in range(K+1):
                    for j in range(K+1-i):
                        nf = f[x][idx][i] + f[y][-1][j] + j*(K-j)
                        if nf > f[x][idx+1][i+j]:
                            f[x][idx+1][i+j] = nf
                            r[x][idx+1][i+j] = j
                idx += 1
    calc_dp(0, -1)
    total_dists = f[0][-1][K]

    chosen = []

    def reconstruct(x, parent, k, chosen):
        idx = len(f[x]) - 1
        for y in reversed(list(E[x])):
            if y != parent:
                assert idx >= 0
                y_k = r[x][idx][k]
                reconstruct(y, x, y_k, chosen)
                k -= y_k
                idx -= 1

        assert idx == 0
        assert 0 <= k and k <= 1
        if k == 1:
            chosen.append(x)

    reconstruct(0, -1, K, chosen)
    return (chosen, total_dists)


def forest_solve(E, C,IndividualToFounderDamMap, RemainingToDoAvailablePerDamLineMap,tree_solver):
    """
    Takes a forest, and a function that solve problem on a tree (e.g. greedy_solve or dp_solve).
    Splits K proportionaly among trees in the forest and solves each of them separately
    using tree_solver. Usually it won't be possible to perfectly split K chosen nodes, so
    number of chosen nodes may be less than K.
    """

    V = list(sorted(E.keys()))
    visited = set()
    chosen = {}
    ChosenSampledPerDamLineMap = {}
    for dam in RemainingToDoAvailablePerDamLineMap:
        ChosenSampledPerDamLineMap[dam] = 0

    total_dist = 0
    visitedFounders = []
    numberOfChosen = 0
    for r in V:
        if r in visited:
            continue

        # a new component
        component = []

        def dfs(x):
            if x not in visited:
                visited.add(x)
                component.append(x)
                for y in E[x]:
                    dfs(y)
        dfs(r)

        if IndividualToFounderDamMap[r] in RemainingToDoAvailablePerDamLineMap:
            K_r = RemainingToDoAvailablePerDamLineMap[IndividualToFounderDamMap[r]]
        else:
            continue

        visitedFounders.append(IndividualToFounderDamMap[r])
        if K_r > 1:
            idxToNode = dict(enumerate(component))
            nodeToIdx = dict([(node, idx) for idx, node in idxToNode.items()])

            # filter out the component edges, but use indexes 0..(N_r-1)
            E_r = [
                list(map(nodeToIdx.get, E[x]))
                for x in component
            ]
            C_r = set(map(nodeToIdx.get, C & set(component)))
            (chosen_r, total_dist_r) = tree_solver(E_r, C_r, K_r)
            total_dist += total_dist_r
            chosen[IndividualToFounderDamMap[r]] = tuple(map(idxToNode.get, chosen_r))
            numberOfChosen += len(tuple(map(idxToNode.get, chosen_r)))
        elif K_r == 1:
            idxToNode = dict(enumerate(component))
            total_dist += 1
            chosen[IndividualToFounderDamMap[r]] = (r,)
            numberOfChosen += 1

    #account for the case of a single individual in line
    if numberOfChosen < sum(RemainingToDoAvailablePerDamLineMap.values()):
        omittedFounders = set(chosen.keys()) ^ set(RemainingToDoAvailablePerDamLineMap.keys())
        for omittedFounder in omittedFounders:
            if (IndividualToFounderDamMap[omittedFounder] == omittedFounder) and (RemainingToDoAvailablePerDamLineMap[omittedFounder] == 1):
                chosen[omittedFounder] = (omittedFounder,)
                numberOfChosen += 1

    return (chosen, total_dist)


def preprocess_input(ParentMap):
    E = defaultdict(list)

    for key, value in ParentMap.items():
        if value == '0':
            continue
        E[value].append(key)
        E[key].append(value)

    return E


"""
Main program reads forest from a given input file.
It outputs a subset of K nodes with maximum total pair-wise distance.
The file format is:
  u_1 v_1
  u_2 v_2
  ..
  u_M v_M
  x1
  x2
  ..
  x_C
Where (u_1, v1_), .., (u_M, v_M) are edges of the forest.
Set {x1, .., x_C}, if included, restricts the choice of solution
nodes to that set. If omitted, any node can be chosen.
"""

def fastSampler(ParentMap,C,IndividualToFounderDamMap,HaplotypedList,RemainingToDoAvailablePerDamLineMap,method,mode):

    solver = None
    if method == 'greedy':
        solver = greedy_solve
    elif method == 'optimal':
        solver = dp_solve
    else:
        assert False, 'Unknown method: {}'.format(method)

    E = preprocess_input(ParentMap)
    C = set(C)
    C.difference_update(set(HaplotypedList))
    start = time.time()
    (chosen, total_dist) = forest_solve(E,C,IndividualToFounderDamMap,RemainingToDoAvailablePerDamLineMap,solver)
    end = time.time()
    elapsed = end - start
    message = f"{method.upper()}:\n\tscore: {total_dist}\n\ttime: {elapsed:8.6f}\n"
    if mode == "gui":
        messagebox.showinfo("Sampling finished", message)
    elif mode == "cl":
        sys.stdout.write("Sampling finished\n")
        sys.stdout.write(message)    

    return chosen
