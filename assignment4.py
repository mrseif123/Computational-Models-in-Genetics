import numpy as np
from matplotlib import pyplot as plt
STOP = 11
N = {1, 3, 6}


def choose_next_successor(h, n, successor_v, values):
    val = max(min(2 * (n - h), 10), 0)
    if h > n:
        values.append(successor_v[val])
    else:
        values.append(h + successor_v[val])


def plot_population(harvests, populations, t):
    x = np.arange(0, t + 1)
    fig = plt.figure()
    fig.set_size_inches(20, 10)
    plt.title("Population")
    plt.ylabel("t")
    plt.plot(x, harvests, '-o')
    plt.plot(x, populations, '-o')
    plt.plot(x, np.subtract(populations, harvests), '-o')
    plt.show()


def update_values(harvests, optimal_h, populations, t):
    populations.append(min(2 * (populations[t - 1] - harvests[t - 1]), 10))
    harvests.append(optimal_h[f"{populations[t]}, {t}"])


def optimize(vals, t):
    current_values = vals
    optimal_harvest_values = {}
    for t in range(t, 0, -1):
        successor_values = {}
        for n in range(STOP):
            values = []
            for h in range(STOP):
                choose_next_successor(h, n, current_values, values)

            successor = max(values)
            optimal_h, optimal_v = values.index(successor), successor
            optimal_harvest_values[f"{n}, {t}"] = optimal_h
            successor_values[n] = optimal_v

        current_values = successor_values
    return optimal_harvest_values


def exercise_1():
    arr = np.arange(0, 10 + 1)
    t = 4
    for n_0 in N:
        optimal_h = optimize(arr, t)
        populations = [n_0]
        harvests = [0]
        for t in range(1, t + 1):
            update_values(harvests, optimal_h, populations, t)
        plot_population(harvests, populations, t)


if __name__ == '__main__':
    exercise_1()
