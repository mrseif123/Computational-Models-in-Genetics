# Computational Models in Genetics and Living Systems
# Home assignment #2, questions 1, 2 & 3.
# By: Seaf Aliyan, mrseif123.

## Question 1, NK Models:-
def process_number(n_as_lst):
    s = ""
    for i in n_as_lst:
        s = s + str(i)
    return s


def generateAllBinaryStrings(n, arr, i, bin_lst):
    if i == n:
        bin_lst.append(process_number(arr[0:n]))
        return

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


def get_neighbours_map(bin_lst, neighbours_map):
    normal_list = get_decimal_representation(bin_lst)
    for binary in bin_lst:
        neighbours_map[binary] = []

    for a in normal_list:
        for b in normal_list:
            if a != b and differAtOneBitPos(a, b):
                tmp_lst = neighbours_map[bin_lst[a]]
                tmp_lst.append(bin_lst[b])
                neighbours_map[bin_lst[a]] = tmp_lst


def simulate_NK(N=14, K=[0, 4, 10]):
    pass


if __name__ == '__main__':
    # Question 1:
    # Setting up environment:-
    N = 4
    helper_list = [None] * N
    bin_lst = []
    neighbours_map = {}
    generateAllBinaryStrings(N, helper_list, 0, bin_lst)
    get_neighbours_map(bin_lst, neighbours_map)
    print(neighbours_map)
