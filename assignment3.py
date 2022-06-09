import numpy as np
import operator as op
import random
from functools import reduce
from matplotlib import pyplot as plt


def Wright_Fisher(N):
    A_s = N * [1]
    a_s = N * [0]
    al = A_s + a_s
    gens = 0
    while True:
        al = np.random.choice(al, size=2 * N, replace=True)
        s = sum(al)
        if s == 0 or s == (2 * N):
            break
        gens += 1

    return gens


def Moran():
    N = 50
    A_s = N * [1]
    a_s = N * [0]
    al = A_s + a_s
    gens = 0
    while True:
        for i in range(len(al)):
            al[i] = np.random.choice(al, replace=True)

        s = sum(al)
        if s == 0 or s == 2 * N:
            break
        gens += 1

    return gens


# print(al)


def Wrigh_Fisher_Model(n):
    GENS_WF = []
    for i in range(500):
        GENS_WF.append(Wright_Fisher(n))
    fig = plt.figure()
    fig.set_size_inches(20, 10)
    plt.hist(x=GENS_WF, rwidth=0.8)
    plt.xlabel("Fixation time")
    plt.ylabel("Count")
    plt.title("Distribution of Fixation Times Wright-Fisher with N={}".format(n))
    plt.show()
    print(np.mean(GENS_WF))


def Moran_Model():
    GENS_M = []
    for i in range(500):
        GENS_M.append(Moran())
    fig = plt.figure()
    fig.set_size_inches(20, 10)
    plt.hist(x=GENS_M, rwidth=0.8)
    plt.xlabel("Fixation time")
    plt.ylabel("Count")
    plt.title("Distribution of Fixation Times Moran")
    plt.show()
    print(np.mean(GENS_M))


from math import comb


def Wright_Fisher_probs(N):
    T = []
    for i in range(N):
        row = []
        for j in range(N):
            p = 1 / (2 * N)
            q = 1 - p
            ncr = comb(2 * N, j)
            T_ij = ncr * (p ** j) * (q ** (2 * N - j))
            row.append(T_ij)
        T.append(row)
    T = np.array(T)
    print("Wright-Fisher model variance is {}".format(T.var()))
    return T


def Moran_probs(N):
    T = []
    for i in range(N):
        row = []
        for j in range(N):
            p = i / 2 * N
            if (j == i + 1):
                T_ij = p * (1 - p)
            elif (j == i - 1):
                T_ij = (1 - p) * p
            elif (j == i):
                T_ij = p ** 2 + (1 - p) ** 2
            else:
                T_ij = 0
            row.append(T_ij)

        T.append(row)
    T = np.array(T)
    T = np.transpose(T)
    print("Moran model variance is {}".format(T.var()))
    return T


if __name__ == '__main__':
    # Question 1:
    # A
    Wrigh_Fisher_Model(50)
    Moran_Model()
    # B
    Wright_Fisher_probs(50)
    Moran_probs(50)
    # C
    N_s = [10, 20, 50, 100]
    for n in N_s:
        Wrigh_Fisher_Model(n)

    # Question 2:
