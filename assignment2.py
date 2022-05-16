# Computational Models in Genetics and Living Systems
# Home assignment #2, questions 1, 2 & 3.
# By: Seaf Aliyan, mrseif123, 211367164.

## Question 1, NK Models:-
def generate_binary_language_helper(N, part, lst):
    if N - 1:
        generate_binary_language_helper(N - 1, part + '0', lst)
        generate_binary_language_helper(N - 1, part + '1', lst)
    else:
        lst.append('1' + part)


def generate_binary_language(lst=[], N=14):
    generate_binary_language_helper(N, '', lst)
    return lst


def simulate_NK(N=14, K=[0, 4, 10]):
    pass


if __name__ == '__main__':
    bin_list = generate_binary_language(lst=[], N=14)
    print(bin_list)