import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats

matplotlib.use('TKAgg')


def get_equi_from_csv():
    f = open('equilibration_time.txt', 'r')
    data = f.read()
    data = data.split(',\n')
    data.pop()
    data = np.array(data)
    data = data.astype(int)
    return data


def gol_histogram(dat):
    hist, bins = np.histogram(dat, bins=35, density=True)

    plt.bar(bins[:-1], hist, width=(bins[1] - bins[0]))
    plt.title('A Histogram of Game of Life Equilibration Time')
    plt.ylabel('Probability Density')
    plt.xlabel('Equilibration Time')
    plt.show()


def get_com():
    with open('com.txt', 'r') as data:
        timestep = []
        xs = []
        ys = []

        for line in data:
            p = line.split()
            timestep.append(p[0])
            xs.append(float(p[1]))
            ys.append(float(p[2]))

    return np.array(timestep), np.array(xs), np.array(ys)


def com_plot(times, xs, ys):
    plt.rcParams['figure.figsize'] = [10, 6]

    times = times.astype(float)

    rs = np.sqrt(xs ** 2  + (ys ** 2))

    # Use linregress to find the gradient and intercept of a fitted line
    res_x = scipy.stats.linregress(times, xs)
    res_y = scipy.stats.linregress(times, ys)
    res_r = scipy.stats.linregress(times, rs)

    # Plotting the x coordinate
    fig, ax = plt.subplots(1, 1)
    ax.scatter(times, xs, label='CoM data')
    ax.plot(times, res_x.slope * times + res_x.intercept, 'r', label='fitted line')

    plt.title('A graph of Centre of Mass x coordinate against time')
    plt.text(0, 21, 'The x-velocity of the glider is ' + str(res_x.slope))
    plt.ylabel('CoM')
    plt.xlabel('Time')
    plt.legend()

    plt.show()

    # plotting the y coordinate
    fig, ax = plt.subplots(1, 1)
    ax.scatter(times, ys, label='CoM data')
    ax.plot(times, res_y.slope * times + res_y.intercept, 'r', label='fitted line')

    plt.title('A graph of Centre of Mass y coordinate against time')
    plt.text(0, 21, 'The y-velocity of the glider is ' + str(res_y.slope))
    plt.ylabel('CoM')
    plt.xlabel('Time')
    plt.legend()

    plt.show()

    # plotting the r vector
    fig, ax = plt.subplots(1, 1)
    ax.scatter(times, rs, label='CoM data')
    ax.plot(times, res_r.slope * times + res_r.intercept, 'r', label='fitted line')

    plt.title('A graph of Centre of Mass y coordinate against time')
    plt.text(0, 30, 'The r-vector velocity of the glider is ' + str(res_r.slope))
    plt.ylabel('CoM')
    plt.xlabel('Time')
    plt.legend()

    plt.show()


def get_infected_data():
    with open('phases.txt', 'r') as data:
        p1_range = []
        p3_range = []
        av_infections = []

        for line in data:
            p = line.split()
            p1_range.append(float(p[0]))
            p3_range.append(float(p[1]))
            av_infections.append(float(p[2]))

    return np.array(p1_range), np.array(p3_range), np.array(av_infections)


def show_phase_diagram(infs):
    p1_range = np.arange(0, 1, 0.05)
    p3_range = np.arange(0, 1, 0.05)

    Z = np.ones((len(p1_range), len(p3_range)))

    counter = 0

    for i in range(len(p1_range)):
        for j in range(len(p3_range)):
            Z[i, j] = infs[counter]
            counter += 1

    Z = np.array(Z)

    plt.contourf(p1_range, p3_range, Z, 1000, cmap='inferno')
    plt.colorbar(label='Fraction Infected')
    plt.title('P1-P3 Phase Diagram')
    plt.xlabel('P1 - Probability of Infection')
    plt.ylabel('P3 - Probability of Becoming Susceptible')
    plt.show()

    plt.contourf(p1_range, p3_range, Z, 2, cmap='inferno')
    plt.colorbar(label='Fraction Infected')
    plt.title('P1-P3 Phase Diagram')
    plt.xlabel('P1 - Probability of Infection')
    plt.ylabel('P3 - Probability of Becoming Susceptible')
    plt.show()


def get_variance_data():
    f = open('variances.txt', 'r')
    data = f.read()
    data = data.split('\n')
    data.pop(0)
    data.pop(-1)

    variances = []
    errors = []
    p1s = []

    for i in data:
        line = i.split(' ')
        p1s.append(float(line[0]))
        variances.append(float(line[1]))
        errors.append(float(line[2]))

    p1s = np.array(p1s)
    variances = np.array(variances)
    errors = np.array(errors)

    return p1s, variances, errors


def show_variance(p1s, variances, errors):
    plt.errorbar(p1s, variances, yerr=errors, ecolor='black', capsize=5)
    plt.title('Plot of Variance vs Infection Probability')
    plt.xlabel('P1 - Probability of Infection')
    plt.ylabel('Variance')
    plt.show()


def get_immune_frac_data():
    f = open('infected_vs_immune_frac.txt', 'r')
    data = f.read()
    data = data.split('\n')
    data.pop(-1)

    immunes = []
    infected = []

    for i in data:
        line = i.split(' ')
        immunes.append(float(line[0]))
        infected.append(float(line[1]))

    return immunes, infected


def plot_immune_frac(immunes, infected):
    plt.plot(immunes, infected)

    plt.title('Graph of Infected Cells vs Immune Cells')
    plt.xlabel('Fraction of Immune Cells')
    plt.ylabel('Fraction of Infected Cells')
    plt.show()


'''
Uncomment function pairs to graph from file.
'''

# timestep, xs, ys = get_com()
# com_plot(timestep, xs, ys)

# data = get_equi_from_csv()
# gol_histogram(data)

# p1s, p3s, infs = get_infected_data()
# show_phase_diagram(infs)

# p1s, variances, errs = get_variance_data()
# show_variance(p1s, variances, errs)

# immune, infected = get_immune_frac_data()
# plot_immune_frac(immune, infected)
