from matplotlib import pyplot as plt

N = range(11)


def find_optimal_single_step(n, v_terminal):  # n = [0.....10]
    i = 0
    v = [-1] * len(n)
    for h in n[::-1]:
        # print(len(n_s), n_s, len(v_terminal), v_terminal, len(H), H, i)
        n_next = min(2 * v_terminal[i], 10) - H[i]  # TODO fix H index
        v_terminal.append(n_next)
        v[h] = h + v_terminal[-1]
        i += 1

    v_opt = max(v)
    h_opt = v.index(v_opt)
    return h_opt, v_opt


def find_optimal_h(v_terminal, time):
    v_next = v_terminal
    H_opt = len(N) * [[-1] * time]
    V_next_new = []
    for t in range(time, 0, -1):
        for n in N:
            h_opt, v_opt = find_optimal_single_step(N, v_next)
            # print(H_opt, n ,t, h_opt)
            H_opt[n][t - 1] = h_opt
            V_next_new.append(v_opt)
        v_terminal.append(V_next_new[-1])
    return H_opt


def find_specific_solution(time, n0):
    n_s = [n0]
    for t in range(1, time + 1):
        x = find_optimal_h(n_s, t)

        # print(n_s[t - 1])
        # print(len(x),x, t)
        H[t - 1] = x[n_s[t - 1] - 1][0]
        n_s.append(min(2 * n_s[t - 1], 10) - H[t])  # TODO fix H index

    return H[:time], n_s[:time]


def plot_pop(pop, t, colors, flag):
    fig = plt.figure()
    fig.set_size_inches(8, 8)
    i = 0
    n = [1, 4, 6]
    for p in pop:
        plt.scatter(x=range(TIME), y=p[:TIME], color=colors[i])
        plt.xlabel("time")
        plt.ylabel("population size ")
        if flag:
            plt.title("Optimal eradication of an\n invasive species with n=10, c=1")
        else:
            plt.title(
                "Population density across different time frames with time=4 and n={}".format(n[i]))
        i = i + 1
        plt.show()


def find_optimal_single_step2(n, v_terminal):  # n = [0.....10]
    i = 0
    v = [-1] * len(n)
    for h in n[::-1]:
        # print(len(n_s), n_s, len(v_terminal), v_terminal, len(H), H, i)
        n_next = v_terminal[i] - R[i]  # TODO fix H index
        v_terminal.append(n_next)
        v[h] = h + v_terminal[-1]
        i += 1

    v_opt = min(v)
    h_opt = v.index(v_opt)
    return h_opt, v_opt


def find_optimal_h2(v_terminal, time):
    v_next = v_terminal
    R_opt = len(N) * [[-1] * time]
    V_next_new = []
    for t in range(time, 0, -1):
        for n in N:
            r_opt, v_opt = find_optimal_single_step2(N, v_next)
            R_opt[n][t - 1] = r_opt
            V_next_new.append(v_opt)
        v_terminal.append(V_next_new[-1])
    return R_opt


def find_specific_solution2(time, n0, c):
    n_s = [n0]
    for t in range(1, time + 1):
        x = find_optimal_h2(n_s, t)
        R[t - 1] = x[-1][0]
        n_s.append(n_s[t] - R[t + 1])  # TODO fix H index

    return R[:time], n_s[:time]


if __name__ == '__main__':
    # Exercise 1:
    TIME = 4

    H = [-1] * len(N)
    res_H, res_n_s1 = find_specific_solution(TIME, 1)
    print(res_H)
    print(res_n_s1)
    print("*************************")

    H = [-1] * len(N)
    res_H, res_n_s2 = find_specific_solution(TIME, 3)
    print(res_H)
    print(res_n_s2)
    print("*************************")

    H = [-1] * len(N)
    res_H, res_n_s3 = find_specific_solution(TIME, 6)
    print(res_H)
    print(res_n_s3)
    print("*************************")

    plot_pop([res_n_s1, res_n_s2, res_n_s3], TIME, ["red", "blue", "green"], False)

    # Exercise 2:
    R = [-1] * len(N)
    TIME = 5
    res_H, res_n_s4 = find_specific_solution2(TIME, 10, 1)
    print(res_H)
    print(res_n_s4)
    plot_pop([res_n_s4], TIME, ["purple"], True)
