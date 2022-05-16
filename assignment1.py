# Computational Models in Genetics and Living Systems
# Home assignment #1, exercises 1, 2 & 3.
# By: Seaf Aliyan, mrseif123, 211367164.
import math
import random
from time import sleep

import matplotlib.pyplot as plt
import numpy as np


def let_one_week_pass(s1, s2, s1_hist, s2_hist, site_1_immigration_rate, site_2_immigration_rate):
    new_site_1 = s1 * (1 - site_1_immigration_rate)
    new_site_1 = new_site_1 + (s2 * site_2_immigration_rate)

    new_site_2 = s2 * (1 - site_2_immigration_rate)
    new_site_2 = new_site_2 + (s1 * site_1_immigration_rate)

    s1 = new_site_1
    s2 = new_site_2

    s1_hist.append(s1)
    s2_hist.append(s2)
    return s1, s2


def loop_n_weeks(site_1, site_2, site_1_weeks_history, site_2_weeks_history, site_1_immigration_rate,
                 site_2_immigration_rate, n=30):
    for i in range(1, n + 1):
        site_1, site_2 = let_one_week_pass(site_1, site_2, site_1_weeks_history, site_2_weeks_history,
                                           site_1_immigration_rate, site_2_immigration_rate)


def print_and_plot_weeks_history(site_1_weeks_history, site_2_weeks_history, title):
    for k in range(len(site_1_weeks_history)):
        print("Week " + str(k) + " = " + str(site_1_weeks_history[k]) + "  " + str(
            site_2_weeks_history[k]))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.title(title)
    plt.ylabel("week")
    plt.xlabel("population count")
    ax1.scatter(site_1_weeks_history, range(len(site_1_weeks_history)), s=10, c='b', marker="s",
                label='site 1')
    ax1.scatter(site_2_weeks_history, range(len(site_1_weeks_history)), s=10, c='r', marker="o",
                label='site 2')
    plt.legend(loc='upper center')
    plt.show()


def simulate_system(site_1_weeks_history, site_2_weeks_history):
    def plot_until_value(n):
        site_1 = []
        site_2 = []
        i = 1
        while i < len(site_1_weeks_history) and \
                (site_1_weeks_history[i] - site_1_weeks_history[i + 1] > n):
            site_1.append(site_1_weeks_history[i] / (site_1_weeks_history[i] + site_2_weeks_history[i]))
            site_2.append(site_2_weeks_history[i] / (site_1_weeks_history[i] + site_2_weeks_history[i]))
            i = i + 1
        print(
            f'Delta={n} | stopped after {i} weeks | n1={site_1_weeks_history[i - 1]} | n2={site_2_weeks_history[i - 1]}'
            f' population structure(s1::s2)={site_1_weeks_history[i] / (site_1_weeks_history[i] + site_2_weeks_history[i])} :: {site_2_weeks_history[i] / (site_1_weeks_history[i] + site_2_weeks_history[i])}')
        plot_overlapping(site_1, site_2,
                         "Population structure raised to the extreme to see it's equilibrium\n"
                         " waiting for difference of less than " + str(n) +
                         " stopped after " + str(i) + " iterations")

    plot_until_value(0.01)
    plot_until_value(0.001)
    plot_until_value(0.000001)
    plot_until_value(0.000000001)
    plot_until_value(0.000000000001)
    plot_until_value(0.000000000000001)


def plot_overlapping(site_1_weeks_history, site_2_weeks_history, title, site_1_weeks_history_matrix=None,
                     site_2_weeks_history_matrix=None, ):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.title(title)
    plt.ylabel("week")

    plt.xlabel("population count")
    if (not site_1_weeks_history_matrix):
        plt.xlabel("population structure")

    ax1.scatter(site_1_weeks_history, range(len(site_1_weeks_history)), s=15, c='b', marker="s",
                label='site 1')
    ax1.scatter(site_2_weeks_history, range(len(site_1_weeks_history)), s=15, c='r', marker="o",
                label='site 2')
    if (site_1_weeks_history_matrix):
        ax1.scatter(site_1_weeks_history_matrix, range(len(site_1_weeks_history)), s=5, c='y', marker="s",
                    label='site 1 matrix')
        ax1.scatter(site_2_weeks_history_matrix, range(len(site_1_weeks_history)), s=5, c='g', marker="o",
                    label='site 2 matrix')
    plt.legend(loc='upper center')
    plt.show()


def loop_one_week_matrix(transition_matrix, n, site_1_weeks_history, site_2_weeks_history):
    res = np.multiply(transition_matrix, n)
    print(res)

    site_1_weeks_history.append(res[0][0] + res[0][1])
    site_2_weeks_history.append(res[1][0] + res[1][1])


def loop_30_weeks_matrix(transition_matrix, n, site_1_weeks_history, site_2_weeks_history):
    for i in range(1, 31):
        powered_matrix = np.linalg.matrix_power(transition_matrix, i)
        loop_one_week_matrix(powered_matrix, n, site_1_weeks_history, site_2_weeks_history)


def ex2_A_B(time, one_plot=False):
    if one_plot:
        r_s = [1.0, 1.1, 1.4, 2, 3, 4]
        results_r = []
        results_a = []

        for i in range(len(r_s)):
            t = 0

            while t < time:
                n_t = r_s[i] ** t
                results_r.append(n_t)
                t = t + 1

            a_s = np.log([1.0, 1.1, 1.4, 2, 3, 4])
            t = 0
            while t < time:
                n_t = np.exp(a_s[i] * t)
                results_a.append(n_t)
                t = t + 1
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.title("Discrete & Continuous linear growth simulation \n multiple graphs on the same plot")
        plt.ylabel("log population count")
        plt.xlabel("time")

        ax1.scatter(results_r, np.log(range(len(results_r))), s=30, c='r', marker="s",
                    label='r discrete')
        ax1.scatter(results_a, np.log(range(len(results_r))), s=5, c='b', marker="s",
                    label='a continuous')

        plt.legend(loc='lower right')
        plt.show()

    else:
        r_s = [1.0, 1.1, 1.4, 2, 3, 4]
        for i in range(len(r_s)):
            results_r = []
            t = 0

            while t < time:
                n_t = r_s[i] ** t
                results_r.append(n_t)
                t = t + 1

            a_s = np.log([1.0, 1.1, 1.4, 2, 3, 4])
            results_a = []
            t = 0
            while t < time:
                n_t = np.exp(a_s[i] * t)
                results_a.append(n_t)
                t = t + 1
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            plt.title("Discrete time simulation with r = " + str(r_s[i]) + " & a = " + str(a_s[i]) +
                      "\n in " + str(time) + " time unit")
            plt.xlabel("time")
            plt.ylabel("population count")
            ax1.scatter(results_r, range(len(results_r)), s=30, c='r', marker="s",
                        label='r discrete')
            ax1.scatter(results_a, range(len(results_a)), s=5, c='b', marker="s",
                        label='a continuous')
            plt.legend(loc='lower right')
            plt.show()


def ex2_C(time, K=100):
    r_s = [1.0, 1.1, 1.4, 2, 3, 4]
    a_s = np.log([1.0, 1.1, 1.4, 2, 3, 4])
    for i in range(len(r_s)):
        r = r_s[i]
        a = a_s[i]
        n0 = 0.1
        Dt = 0.01
        T = time
        time_steps = math.floor(T / Dt)
        n_a = np.zeros(time_steps + 1)
        n_a[0] = n0
        n_r = np.zeros(time_steps + 1)
        n_r[0] = n0

        for t in range(time_steps):
            n_r[t + 1] = r * n_r[t] * (1 - n_r[t] / K)
            n_a[t + 1] = n_a[t] + Dt * a * (n_a[t] - n_a[t] ** 2 / K)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.title("Logistic dynamics simulation with r = " + str(r_s[i]) + " & a = " + str(a_s[i]) +
                  "\n in " + str(range(len(n_r))) + " time steps")
        plt.ylabel("population count")
        plt.xlabel("time")
        ax1.scatter(range(len(n_r)), n_r, s=7, c='y', marker="s",
                    label='r discrete')
        ax1.scatter(range(len(n_a)), n_a, s=1, c='g', marker="s",
                    label='a continuous')

        plt.legend(loc='upper center')
        plt.show()


def ex3_a(T, n, n0, beta, delta, gamma):
    t = 0
    Dt = 0.01
    times = []
    while t < T:
        t = t + Dt
        for i in range(n):
            rand = random.uniform(0, 1)
            if rand > (1 - beta * Dt):
                n0 = n0 + 1
            elif rand < (1 - delta * Dt):
                n0 = n0 - 1
        for i in range(n * (n - 1)):
            rand = random.uniform(0, 1)
            if rand > (1 - gamma * Dt):
                n0 = n0 - 1
        times.append([t, n0])
        print(t, n0)


def get_ext_time(T, beta, delta, gamma, n0=10):
    res = []
    for i in range(50):
        res.append(gillespie_algo(T, beta, delta, gamma, n0=10))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.title("Logistic dynamics simulation with ")
    plt.ylabel("population size")
    plt.xlabel("time")
    ax1.scatter(res, range(len(res)), s=7, c='y', marker="s",
                label='r discrete')
    plt.legend(loc='upper center')
    plt.show()


def gillespie_algo(T, beta, delta, gamma, n0=10):
    t = 0
    while n0 != 0:
        p_birth = beta
        p_death = delta
        p_comp = gamma * n0 * (n0 - 1)
        p_total = p_birth + p_death + p_comp

        r1 = random.uniform(0, 1)
        sleep(0.001)
        r2 = random.uniform(0, 1)

        if r1 < p_birth / p_total:
            n0 = n0 + 1
        elif r1 < (p_birth + p_death) / p_total:
            n0 = n0 - 1
        else:
            n0 = n0 - 1

        t = t + (1 / p_total) * (np.log(1 / r2))
    print(t, n0)
    return t


if __name__ == '__main__':
    # # ex1 :
    # site_1_weeks_history = []
    # site_2_weeks_history = []
    #
    # site_1 = 100
    # site_2 = 100
    #
    # site_1_weeks_history.append(site_1)
    # site_2_weeks_history.append(site_2)
    #
    # site_1_immigration_rate = 0.18
    # site_2_immigration_rate = 0.12
    #
    # # ex1 : A
    # loop_n_weeks(site_1, site_2, site_1_weeks_history, site_2_weeks_history, site_1_immigration_rate,
    #              site_2_immigration_rate, 30)
    #
    # print_and_plot_weeks_history(site_1_weeks_history, site_2_weeks_history,
    #                              "30 first weeks population of sites 1 & 2")
    #
    # # ex1: B
    # site_1 = 100
    # site_2 = 100
    # site_1_weeks_history_matrix = [site_1]
    # site_2_weeks_history_matrix = [site_2]
    # n_0 = np.array([[site_1],
    #                 [site_2]])
    # m_transition_matrix = np.array([[1 - site_1_immigration_rate, site_2_immigration_rate],
    #                                 [site_1_immigration_rate, 1 - site_2_immigration_rate]])
    #
    # loop_30_weeks_matrix(m_transition_matrix, n_0, site_1_weeks_history_matrix,
    #                      site_2_weeks_history_matrix)
    # print_and_plot_weeks_history(site_1_weeks_history_matrix, site_2_weeks_history_matrix,
    #                              "30 first weeks population of sites 1 & 2,"
    #                              " using transition matrix")
    #
    # plot_overlapping(site_1_weeks_history, site_2_weeks_history, "30 first weeks population of sites 1 & 2,\n"
    #                                                              " using transition matrix and equation calulation"
    #                                                              " overlapped", site_1_weeks_history_matrix,
    #                  site_2_weeks_history_matrix)
    #
    # # ex1: C
    # site_1_weeks_history = []
    # site_2_weeks_history = []
    #
    # loop_n_weeks(site_1, site_2, site_1_weeks_history, site_2_weeks_history, site_1_immigration_rate,
    #              site_2_immigration_rate, 5000)
    #
    # simulate_system(site_1_weeks_history, site_2_weeks_history)
    #
    # # ex1: D (Bonus)
    # w, v = np.linalg.eig(m_transition_matrix)
    # max_eig_v = max(w)
    # print(max_eig_v)
    #
    # # ex2: A & B
    # ex2_A_B(time=50, one_plot=True)
    #
    # # ex2: C
    ex2_C(time=100, K=100)
    #
    # # ex2: D (Bonus)

    # ex3: A
    gillespie_algo(50, 2, 1, 0.1, 10)
    get_ext_time(50, 2, 1, 0.1, 10)
    ex3_a(10, 10, 10, 2, 1, 0.1)
