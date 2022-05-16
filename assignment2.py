# Computational Models in Genetics and Living Systems
# Home assignment #2, questions 1, 2 & 3.
# By: Seaf Aliyan, mrseif123.

## Question 1, NK Models:-
import random

from matplotlib import pyplot as plt

MARKERS_SHAPES = ["v", "^", "8", "s", "p", "*", "h",
                  "+", "X", "d", "D", "<", ">", "1",
                  "2", "3", "4", "."]


def process_number(n_as_lst):
    s = ""
    for i in n_as_lst:
        s = s + str(i)
    return s


def generateAllBinaryStrings(n, arr, i, bin_lst):
    if i == n:
        bin_lst.append(process_number(arr[0:n]))
        return bin_lst

    arr[i] = 0
    generateAllBinaryStrings(n, arr, i + 1, bin_lst)

    arr[i] = 1
    generateAllBinaryStrings(n, arr, i + 1, bin_lst)


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


def update_map(n_m):
    seen_points = set()
    new_map = {}
    for k in n_m.keys():
        if k not in seen_points:
            seen_points.add(k)
            for v in n_m[k]:
                seen_points.add(v)
            new_map[k] = n_m[k]
    n_m = new_map


def create_shapes_map(n_m):
    shapes_map = {}
    seen_points = set()
    counter = 0
    for k in n_m.keys():
        shapes_map[k] = MARKERS_SHAPES[counter]
        seen_points.add(k)
        for v in n_m[k]:
            if v not in seen_points:
                counter += 1
                seen_points.add(v)
                shapes_map[v] = MARKERS_SHAPES[counter]
            if v in seen_points:
                counter-=1
    return shapes_map


def plot_fitness(fitness_map, N, neighbours_map):
    update_map(neighbours_map)
    # shapes_map = create_shapes_map(neighbours_map)
    keys = fitness_map.keys()
    vals = [fitness_map[i] for i in keys]
    plt.scatter(keys, vals)
    plt.plot(keys, vals)
    plt.xlabel("")
    plt.ylabel("fitness")
    plt.title("Fitness landscape\n with N=" + str(N) + ", K=" + str(K))
    plt.show()
    pass


if __name__ == '__main__':
    # Question 1:
    # Setting up environment:-
    N = 5
    K = 4
    helper_list = [None] * N
    bin_lst = []
    generateAllBinaryStrings(N, helper_list, 0, bin_lst)
    neighbours_map = get_neighbours_map(bin_lst)
    fitness_map = get_local_fitness_lst(N, K, bin_lst)
    plot_fitness(fitness_map, N, neighbours_map)
