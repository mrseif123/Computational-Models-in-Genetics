# Computational Models in Genetics and Living Systems
# Home assignment #2, questions 1, 2 & 3.
# By: Seaf Aliyan, mrseif123.

# Question 1, NK Models:-
import random

import numpy as np
from matplotlib import pyplot as plt
from statsmodels.graphics import tsaplots

MARKERS_SHAPES = [".", "v", "+", "*", "^", "8", "s", "p", "h",
                  "X", "d", "D", "<", ">", "1", "2", "3", "4"]


def process_number(n_as_lst):
    s = ""
    for i in n_as_lst:
        s = s + str(i)
    return s


def generateAllBinaryStrings(num, arr, i, bin_lst):
    if i == num:
        bin_lst.append(process_number(arr[0:num]))
        return bin_lst

    arr[i] = 0
    generateAllBinaryStrings(num, arr, i + 1, bin_lst)

    arr[i] = 1
    generateAllBinaryStrings(num, arr, i + 1, bin_lst)


def isPowerOfTwo(x):
    # First x in the below expression is
    # for the case when x is 0
    return x and (not (x & (x - 1)))


# function to check whether the two numbers
# differ at one bit position only
def differAtOneBitPos(a, b):
    return isPowerOfTwo(a ^ b)


def dec_rep(n):
    return int(n, 2)


def get_decimal_representation(lst):
    decimal_list = list(map(dec_rep, lst))
    return decimal_list


def get_neighbours_map(bin_lst):
    neighbours_map = {}
    normal_list = get_decimal_representation(bin_lst)
    for binary in bin_lst:
        neighbours_map[binary] = []

    for a in normal_list:
        for b in normal_list:
            if a != b and differAtOneBitPos(a, b):
                tmp_lst = neighbours_map[bin_lst[a]]
                tmp_lst.append(bin_lst[b])
                neighbours_map[bin_lst[a]] = tmp_lst

    return neighbours_map


def generate_fitness_uniform(K):
    tmp = [None] * (K + 1)
    binary_of_length_k = []
    fitness = {}
    generateAllBinaryStrings(K + 1, tmp, 0, binary_of_length_k)
    for i in binary_of_length_k:
        fitness[i] = random.uniform(0, 1)
    return fitness, binary_of_length_k


def sum_fitness_value(fi, K, N, fi_s):
    val = 0
    tmp_fi = fi + fi + fi
    for i in range(N):
        index = (i + K + 1)
        val = val + fi_s[tmp_fi[i: index]]
    return val


def get_local_fitness_lst(N, K, bin_lst):
    fi_s, bi_s = generate_fitness_uniform(K)
    fitness_map = {}
    for b in bin_lst:
        fitness_map[b] = 0

    for fi in fitness_map.keys():
        fitness_map[fi] = float(sum_fitness_value(fi, K, N, fi_s) / N)

    return fitness_map


def create_shapes_map(n_m):
    shapes_map = {}
    for k in n_m.keys():
        shapes_map[k] = MARKERS_SHAPES[k.count('1')]
    return shapes_map


def plot_fitness(fitness_map, N, neighbours_mapx):
    shapes_map = create_shapes_map(neighbours_map)
    keys = fitness_map.keys()
    vals = [fitness_map[i] for i in keys]

    key_lst = list(keys)
    for i in range(len(key_lst)):
        plt.scatter(key_lst[i], vals[i])
    for k in neighbours_mapx.keys():
        my_vals = [fitness_map[v] for v in neighbours_mapx[k]]
        plt.plot(neighbours_mapx[k], my_vals, marker=shapes_map[key_lst[i]])
    plt.xlabel("")
    plt.ylabel("fitness")
    plt.title("Fitness landscape\n with N=" + str(N) + ", K=" + str(K))
    plt.xticks([])
    plt.show()
    pass


def get_trajectory_of_length(n, n_m, strt_point):
    ln = pow(2, n)
    vec = [strt_point]
    curr = strt_point
    for i in range(ln - 1):
        random_neighbours = n_m[curr]
        chosen = random.choice(random_neighbours)
        vec.append(chosen)
    return vec


def calc_fitness(lst, f_map):
    vec = [f_map[x] for x in lst]
    return vec


def autocorr(x):
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


def plot_non_decrasing_trajectories(b_lst, n_m, f_m, k, n):
    t_l = get_trajectory_lengths(b_lst, f_m, n_m)
    fig = plt.figure(figsize=(10, 5))
    plt.bar(t_l.keys(), t_l.values())
    plt.xlabel("starting point")
    plt.ylabel("path length")
    plt.title("Distrubution of longest non-decreasing\n fitness trajectories with N={} & K={}".format(N,K))
    plt.show()


def get_trajectory_lengths(b_lst, f_m, n_m):
    trajectory_lenght_map = {}
    absloute_max = max([val for val in f_m.values()])
    for point in b_lst:
        trajectory_lenght_map[point] = 0
        current_max = point
        current_max_val = f_m[current_max]
        while True:
            neighbours_fitnesses = [f_m[current_max]]
            for neighbour in n_m[current_max]:
                if f_m[neighbour] > current_max_val:
                    current_max_val = f_m[neighbour]
                    current_max = neighbour
                    trajectory_lenght_map[point] += 1
            if current_max_val != absloute_max:
                trajectory_lenght_map[point] += 1
            else:
                break
    return trajectory_lenght_map


def autocorrelation_flow():
    start_point = random.choice(bin_lst)
    trajectory = get_trajectory_of_length(N, neighbours_map, start_point)
    trajectory_as_fitness = np.array(calc_fitness(trajectory, fitness_map))
    auto_corr = autocorr(trajectory_as_fitness)
    fig = tsaplots.plot_acf(auto_corr, lags=N)
    plt.title("Autocorrelation with " + str(N) + " lags")
    plt.ylabel("autocorrelation")
    plt.show()
    print(auto_corr)


def local_maximums_flow():
    num_of_local_maximums = get_local_maximums_number(neighbours_map, fitness_map)
    print(num_of_local_maximums)


def longest_trajectories_flow():
    plot_non_decrasing_trajectories(bin_lst, neighbours_map, fitness_map, K, N)


if __name__ == '__main__':
    # Question 1:
    for n, k in [(14, 0), (14, 4), (14, 10)]:
        # Setting up environment:-
        N = n
        K = k
        N = 3
        K = 2
        helper_list = [None] * N
        bin_lst = []
        generateAllBinaryStrings(N, helper_list, 0, bin_lst)
        neighbours_map = get_neighbours_map(bin_lst)
        fitness_map = get_local_fitness_lst(N, K, bin_lst)
        plot_fitness(fitness_map, N, neighbours_map) # TODO (1) Check graphs -too crowded for required N=14-

        # Part i:-
        autocorrelation_flow() # TODO (1) Answer questions.

        # Part ii:-
        local_maximums_flow()

        # Part iii:-
        longest_trajectories_flow() # TODO (1) fix infinite loop (2) Answer questions.
        break
