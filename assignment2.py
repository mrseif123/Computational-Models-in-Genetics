# Computational Models in Genetics and Living Systems
# Home assignment #2, questions 1, 2 & 3.
# By: Seaf Aliyan, mrseif123.

# Question 1, NK Models:-
import random

import numpy as np
from matplotlib import pyplot as plt

MARKERS_SHAPES = [".", "v", "+", "*", "^", "8", "s", "p", "h",
                  "X", "d", "D", "<", ">", "1", "2", "3", "4"]

NUM_OF_MAXIMUMS_Q1 = []
NUM_OF_MAXIMUMS_Q2 = -1

MAX_PATH = -1

AUTOCORRELATION_Q1 = []
CORRELATION_Q2 = -1

TRAJECTORIES_DISTRIBUTION_Q1 = -1
TRAJECTORIES_DISTRIBUTION_Q2 = -1


def process_number(n_as_lst):
    s = ""
    for i in n_as_lst:
        s = s + str(i)
    return s


def generateAllBinaryStrings(num, arr, i, bin_list):
    if i == num:
        bin_list.append(process_number(arr[0:num]))
        return bin_list

    arr[i] = 0
    generateAllBinaryStrings(num, arr, i + 1, bin_list)

    arr[i] = 1
    generateAllBinaryStrings(num, arr, i + 1, bin_list)


def isPowerOfTwo(x):
    # First x in the below expression is for the case when x is 0
    return x and (not (x & (x - 1)))


# function to check whether the two numbers differ at one bit position only
def differAtOneBitPos(a, b):
    return isPowerOfTwo(a ^ b)


def dec_rep(n):
    return int(n, 2)


def get_decimal_representation(lst):
    decimal_list = list(map(dec_rep, lst))
    return decimal_list


def get_neighbours_map(bin_list):
    neighbours_m = {}
    normal_list = get_decimal_representation(bin_list)
    for binary in bin_list:
        neighbours_m[binary] = []

    for a in normal_list:
        for b in normal_list:
            if a != b and differAtOneBitPos(a, b):
                tmp_lst = neighbours_m[bin_list[a]]
                tmp_lst.append(bin_list[b])
                neighbours_m[bin_list[a]] = tmp_lst

    return neighbours_m


def generate_fitness_uniform(k):
    tmp = [None] * (k + 1)
    binary_of_length_k = []
    fitness = {}
    generateAllBinaryStrings(k + 1, tmp, 0, binary_of_length_k)
    for i in binary_of_length_k:
        fitness[i] = random.uniform(0, 1)
    return fitness, binary_of_length_k


def sum_fitness_value(fi, k, n, fi_s):
    val = 0
    tmp_fi = fi + fi + fi
    for i in range(n):
        index = (i + k + 1)
        val = val + fi_s[tmp_fi[i: index]]
    return val


def get_local_fitness_lst(n, k, bin_list):
    fi_s, bi_s = generate_fitness_uniform(k)
    f_m = {}
    for b in bin_list:
        f_m[b] = 0

    for fi in f_m.keys():
        f_m[fi] = float(sum_fitness_value(fi, k, n, fi_s) / n)

    return f_m


def get_local_fitness_lst_total_random(n, bin_list):
    f_m = {}
    for b in bin_list:
        f_m[b] = random.randint(0, 1)

    return f_m


def get_trajectory_of_length(n, n_m, s_p):
    ln = pow(2, n)
    vec = [s_p]
    curr = s_p
    for i in range(ln - 1):
        random_neighbours = n_m[curr]
        chosen = random.choice(random_neighbours)
        vec.append(chosen)
    return vec


def calc_fitness(lst, f_map):
    vec = [f_map[x] for x in lst]
    return vec

# Choose the 3rd autocorrelation equation (normalized to 1 at zero):-
def autocorr(x):
    data = x
    lags = range(11)
    acorr = len(lags) * [0]

    # Mean
    mean = sum(data) / len(data)

    # Variance
    var = sum([(x - mean) ** 2 for x in data]) / len(data)

    # Normalized data
    ndata = [x - mean for x in data]

    # Go through lag components one-by-one
    for l in lags:
        c = 1  # Self correlation

        if (l > 0):
            tmp = [ndata[l:][i] * ndata[:-l][i]
                   for i in range(len(data) - l)]

            c = sum(tmp) / len(data) / var

        acorr[l] = c
    return acorr

def corr(x):
    result = np.correlate(x, x, mode='full')
    return result[int(len(result) / 2):]


def get_local_maximums_number(n_m, f_m):
    counter = 0
    for k in n_m.keys():
        flag = True
        for neighbour in n_m[k]:
            if f_m[k] < f_m[neighbour]:
                flag = False
                break
        if flag:
            counter += 1
    return counter


def plot_non_decreasing_trajectories(b_lst, n_m, f_m, k, n):
    t_l, max_path = get_trajectory_lengths(b_lst, f_m, n_m)
    val_as_list = list(t_l.values())
    plt.hist(x=t_l.values(), bins=np.arange(min(val_as_list) + 0.25, max(val_as_list) + 1) + 0.2,
             rwidth=0.6, )

    plt.xlabel("path length")
    plt.ylabel("count")
    plt.xticks(range(1, max_path + 1))
    plt.title(
        "Distribution of longest non-decreasing\n fitness trajectories with N={} & K={}"
        " & number of paths={}".format(n, k, len(t_l.keys())))
    plt.show()


def get_trajectory_lengths(b_lst, f_m, n_m):
    global MAX_PATH

    trajectory_length_map = {}
    absolute_max = max([val for val in f_m.values()])
    for point in b_lst:
        trajectory_length_map[point] = 0
        current_max = point
        current_max_val = f_m[current_max]
        while True:
            for neighbour in n_m[current_max]:
                if f_m[neighbour] > current_max_val:
                    current_max_val = f_m[neighbour]
                    current_max = neighbour
                    trajectory_length_map[point] += 1

            someone_bigger = False
            for neighbour in n_m[current_max]:
                if f_m[neighbour] > current_max_val:
                    someone_bigger = True
            if current_max_val != absolute_max:
                trajectory_length_map[point] += 1
            if not someone_bigger:
                break
    MAX_PATH = max(trajectory_length_map.values())
    return trajectory_length_map, MAX_PATH


def relation_flow(n, k, Q1):
    global AUTOCORRELATION_Q1, CORRELATION_Q2
    start_point = random.choice(bin_lst)
    trajectory = get_trajectory_of_length(N, neighbours_map, start_point)
    trajectory_as_fitness = np.array(calc_fitness(trajectory, fitness_map))

    if Q1:
        relation = autocorr(trajectory_as_fitness)
    else:
        relation = corr(trajectory_as_fitness)

    if Q1:
        AUTOCORRELATION_Q1.append(relation)
        plt.title("Autocorrelation N={} K={}".format(n,k))
        plt.ylabel("autocorrelation")
        rel = relation[::-1] + relation[1::]
        plt.plot( range(-len(rel)//2 + 1, len(rel)//2 + 1),rel)
    else:
        CORRELATION_Q2 = relation
        plt.title("Correlation with N={} K={}".format(n,k))
        plt.ylabel("correlation")

    plt.show()



def local_maximums_flow(Q1):
    global NUM_OF_MAXIMUMS_Q1, NUM_OF_MAXIMUMS_Q2

    if Q1:
        NUM_OF_MAXIMUMS_Q1.append(get_local_maximums_number(neighbours_map, fitness_map))
    else:
        NUM_OF_MAXIMUMS_Q2 = get_local_maximums_number(neighbours_map, fitness_map)


def longest_trajectories_flow():
    plot_non_decreasing_trajectories(bin_lst, neighbours_map, fitness_map, K, N)


if __name__ == '__main__':
    # Question 1:
    print("Question 1 answers:-")

    m = 0
    for N, K in [(14, 0), (14, 4), (14, 10)]:
        # Setting up environment:-
        print("\tStarting for N={}, K={}".format(N, K))

        helper_list = [None] * N
        bin_lst = []
        generateAllBinaryStrings(N, helper_list, 0, bin_lst)
        neighbours_map = get_neighbours_map(bin_lst)
        fitness_map = get_local_fitness_lst(N, K, bin_lst)

        # Part i:-
        relation_flow(N, K, Q1=True)  # TODO (1) Answer questions.
        print("\tAutocorrelation for N={}, K={}  is: {}".format(N, K, AUTOCORRELATION_Q1[m])) # TODO FIX

        # Part ii:-
        local_maximums_flow(Q1=True)
        print("\tNumber of local maximums for N={}, K={}  is: {}".format(N, K, NUM_OF_MAXIMUMS_Q1[m]))

        # Part iii:-
        longest_trajectories_flow()  # TODO (1) fix infinite loop (2) Answer questions.
        print("\tMax path length with N={}, K={} is: {}".format(N, K, MAX_PATH))
        m += 1

    print("######################################################################\n"
          "\n######################################################################")

    # Question 2: # TODO (2.i) Compare them to NK landscape.
    print("Question 2 answers:-\n (this is a private-case for NK model where K=N-1)")

    N = 10
    K = N - 1
    helper_list = [None] * N
    bin_lst = []
    generateAllBinaryStrings(N, helper_list, 0, bin_lst)
    neighbours_map = get_neighbours_map(bin_lst)
    fitness_map = get_local_fitness_lst(N, K, bin_lst)

    # Part i:-
    relation_flow( N, K, Q1=False)
    print("\tCorrelation for N={}, K={} is: {}".format(N, K, CORRELATION_Q2)) # TODO FIX

    # Part ii:-
    local_maximums_flow(Q1=False)
    print("\tNumber of local maximums for N={}, K={} is: {}".format(N, K, NUM_OF_MAXIMUMS_Q2))

    # Part iii:-
    longest_trajectories_flow()
    print("\tMax path length with N={}, K={} is: {}".format(N, K, MAX_PATH))

    # Comparing with NK model:- // TODO Add comparesion in the writeup.
