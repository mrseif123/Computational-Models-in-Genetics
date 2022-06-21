import random
from math import comb
import numpy as np
from matplotlib import pyplot as plt


def Wright_Fisher(N):
    x = N * [1]
    y = N * [0]
    x_y = x + y
    result = 0
    while True:
        x_y = np.random.choice(x_y, size=2 * N, replace=True)
        current_sum = sum(x_y)
        if current_sum == 0 or current_sum == (2 * N):
            break
        result = result + 1
    return result


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


def Wrigh_Fisher_Model(n):
    results = []
    for i in range(500):
        results.append(Wright_Fisher(n))
    fig = plt.figure()
    fig.set_size_inches(20, 10)
    plt.hist(x=results, rwidth=0.8)
    plt.xlabel("Fixation Time")
    plt.ylabel("Count")
    plt.title("Distribution Of Fixation Times Wright Fisher Model with N={}".format(n))
    plt.show()
    print("The Mean of Wright Fisher Model is: {} ".format(np.mean(results)))


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
            if j == i + 1:
                T_ij = p * (1 - p)
            elif j == i - 1:
                T_ij = (1 - p) * p
            elif j == i:
                T_ij = p ** 2 + (1 - p) ** 2
            else:
                T_ij = 0
            row.append(T_ij)

        T.append(row)
    T = np.array(T)
    T = np.transpose(T)
    print("Moran model variance is {}".format(T.var()))
    return T


def Drift_Selection_Wright_Fisher(pops, repeats, w_A, w_a, p):
    for pop in pops:
        m = 0
        for loop in repeats:
            A_s = int(pop / 2) * [1]
            a_s = int(pop / 2) * [0]
            al = A_s + a_s

            while True:
                al_set = set(al)
                al = np.random.choice(al, size=pop, replace=True)
                num_of_zeros = sum([1 for x in al if x == 0])
                num_of_ones = sum([1 for x in al if x == 1])
                w_bar = ((1.1 * num_of_ones) + (1 * num_of_zeros)) / (num_of_ones + num_of_zeros)

                s = sum(al)
                if s == 0 or s == pop:
                    break

            m = m + 1

        print("After {} repeats, mean time to fixation is {}".format(repeats, m))

    return


def ALLESLES_AFTER_N_GENERATIONS(gens, p, n):
    a, b, c, lst = Wright_Fisher_MUTATED(n, p, gens)

    H_0 = 1 / 2
    H_s = [H_0]
    for i in range(gens):
        H_s.append(H_s[i] + (1 - H_s[i]) * (1 - (1 - p) ** 2))

    fig = plt.figure()
    fig.set_size_inches(8, 8)
    plt.scatter(x=range(gens + 1), y=H_s)
    plt.xlabel("Time")
    plt.ylabel("Heterozygotes")
    plt.title("Heterozygotes of Allele at time T, ran for {} generations".format(gens))
    plt.show()

    fig = plt.figure()
    fig.set_size_inches(8, 8)
    plt.bar(x=["A", "a"], height=[b, c], color=["green", "blue"])
    plt.xlabel("Allele Type")
    plt.ylabel("Frequency")
    plt.title("Frequency of Allele At The End of Simulation\n With Mutation Rate of {}"
              " and Population Size {}".format(p, n))
    plt.show()


def Wright_Fisher_MUTATED(n, p, run_gens):
    A_s = n * [1]
    a_s = n * [0]
    all_ale = A_s + a_s
    gens = 0
    r = [0, p * 100]
    counter = 0
    new_ls = []

    while counter < run_gens:
        for a in all_ale:
            rnd_int = random.uniform(r[0], r[1])
            # print(rnd_int)
            if rnd_int < p:
                if a == 1:
                    new_ls.append(0)
                else:
                    new_ls.append(1)
            else:
                new_ls.append(a)

        s = sum(new_ls)
        if s == 0 or s == (2 * n):
            break
        all_ale = new_ls
        new_ls = []
        gens += 1
        counter += 1

    A__s = all_ale.count(1)
    a__s = all_ale.count(0)
    print("Stopped after {} generations \t A number is {} \t a number is {}".format(gens, A__s, a__s))
    return gens, A__s, a__s, all_ale


def calculate_w_bar(ls):
    number_A = ls.count(1)
    number_a = ls.count(0)

    mean = (1.1 * number_A + number_a * 1) / (number_a + number_A)
    return mean


def Wright_Fisher_helper(n, lst, prob):
    al = lst
    gens = 0
    while 100:
        p_s = []
        for i in al:
            if i == 1:
                p_s.append(prob / list(al).count(1))
            else:
                p_s.append((1 - prob) / list(al).count(0))
        al = np.random.choice(al, size=2 * n, replace=True)
        s = sum(al)
        if s == 0 or s == (2 * n):
            break
        gens += 1

    return al


def get_lst(base, w_A, w_a, n):
    prob = []
    base = list(base)
    count_a = base.count(0)
    count_A = base.count(1)

    for i in base:
        if i == 1:
            prob.append((w_A / (w_a + w_A)) / count_A)
        else:
            prob.append((w_a / (w_a + w_A)) / count_a)

    lst = np.random.choice(base, size=2 * n, replace=True, p=prob)

    return lst


def Wright_Fisher_Selection(n, w_A, w_a, prob):
    A_s = (int(n * prob)) * [1]
    a_s = (n - len(A_s)) * [0]
    base = A_s + a_s

    while True:
        lst = get_lst(base, w_A, w_a, n)
        s = sum(lst)
        print(s)
        if s == 0 or s == 2 * n:
            break
        base = lst

    print(s, len(lst), lst)
    if sum(lst) == len(lst):
        print("List fixed to 1")
    elif sum(lst) == 0:
        print("List fixed to 0")
    else:
        print("List didn't fix")
    return lst


if __name__ == '__main__':
    # Question 1:
    # A
    Wrigh_Fisher_Model(50)
    Moran_Model()
    # B
    Wright_Fisher_probs(50)
    Moran_probs(50)
    # C
    N_s = [10, 30, 50, 100]
    for n in N_s:
        Wrigh_Fisher_Model(n)


    # Question 2:
    N = 1000
    GENs = 5000
    MUTATION_RATES = [0.1 / N, 1 / N, 10 / N]
    # i, ii, iii (Different Allele types):
    for p in MUTATION_RATES:
        ALLESLES_AFTER_N_GENERATIONS(GENs, p, N)


    # Question 3:
    REPEATS = 5
    P = 0.02
    W_a = 1
    W_A = 1.1
    POPULATION_SIZES = [51, 1000, 10000]
    for pop_size in POPULATION_SIZES:
        Wright_Fisher_Selection(pop_size, W_A, W_a, P)
