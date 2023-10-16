import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def read_from_file(filename):
    with open(filename, 'r') as data:
        first_col = []
        second_col = []

        for line in data:
            row = line.split(' ')
            first_col.append(float(row[0]))
            second_col.append(float(row[1]))

    return np.array(first_col), np.array(second_col)


def Exponential(x, A, B):
    y = -A * B**x
    return y


def plot_free_energy():
    timestep, energies = read_from_file('free_energy.txt')

    plt.plot(timestep, energies)
    plt.title('Free energy against time')
    plt.xlabel('time')
    plt.ylabel('Free energy')
    plt.show()


def plot_poisson():
    distance, potential = read_from_file('poisson_potential.txt')
    distance_e, electric = read_from_file('electric_field_poisson.txt')

    parameters, covariance = curve_fit(Exponential, distance, potential)
    opt_distance = parameters[0]
    opt_potential = parameters[1]

    fit_potential = Exponential(distance, opt_distance, opt_potential)

    plt.plot(distance, potential, marker='x', label='data')
    plt.plot(distance, fit_potential, '-', label='fit')
    plt.title('Potential strength against distance from charge')
    plt.xlabel('Distance')
    plt.ylabel('Field Strength')
    plt.legend()
    plt.show()

    e_parameters, e_covariance = curve_fit(Exponential, distance, electric)
    opt_e_distance = e_parameters[0]
    opt_electric = parameters[1]
    fit_electric = Exponential(distance, opt_e_distance, opt_electric)

    plt.plot(distance_e, electric, marker='x', label='data')
    plt.plot(distance_e, fit_electric, '-', label='fit')
    plt.title('Electric Field strength against distance from charge')
    plt.xlabel('Distance')
    plt.ylabel('E Field Strength')
    plt.legend()
    plt.show()


def plot_gs():
    distance, potential = read_from_file('gs_potential.txt')
    distance_e, electric = read_from_file('electric_field_gs.txt')

    parameters, covariance = curve_fit(Exponential, distance, potential)
    opt_distance = parameters[0]
    opt_potential = parameters[1]

    fit_potential = Exponential(distance, opt_distance, opt_potential)

    plt.plot(distance, potential, marker='x', label='data')
    plt.plot(distance, fit_potential, '-', label='fit')
    plt.title('Potential strength against distance from charge using \n the Gauss-Seidel Algorithm')
    plt.xlabel('Distance')
    plt.ylabel('Field Strength')
    plt.legend()
    plt.show()

    e_parameters, e_covariance = curve_fit(Exponential, distance, electric)
    opt_e_distance = e_parameters[0]
    opt_electric = parameters[1]
    fit_electric = Exponential(distance, opt_e_distance, opt_electric)

    plt.plot(distance_e, electric, marker='x', label='data')
    plt.plot(distance_e, fit_electric, '-', label='fit')
    plt.title('Electric Field strength against distance from charge using \n the Gauss-Seidel Algorithm')
    plt.xlabel('Distance')
    plt.ylabel('E Field Strength')
    plt.legend()
    plt.show()


def plot_magnetic():
    distance, potential = read_from_file('magnetic_potential.txt')
    distance_b, magnetic = read_from_file('magnetic_field.txt')

    parameters, covariance = curve_fit(Exponential, distance, potential)
    opt_distance = parameters[0]
    opt_potential = parameters[1]

    fit_potential = Exponential(distance, opt_distance, opt_potential)

    plt.plot(distance, potential, marker='x', label='data')
    plt.plot(distance, fit_potential, '-', label='fit')
    plt.title('Potential strength against distance from wire')
    plt.xlabel('Distance')
    plt.ylabel('Field Strength')
    plt.legend()
    plt.show()

    b_parameters, b_covariance = curve_fit(Exponential, distance, magnetic)
    opt_b_distance = b_parameters[0]
    opt_magnetic = b_parameters[1]
    fit_magnetic = Exponential(distance, opt_b_distance, opt_magnetic)

    plt.plot(distance_b, magnetic, marker='x', label='data')
    plt.plot(distance_b, fit_magnetic, '-', label='fit')
    plt.title('Magnetic Field strength against distance from wire')
    plt.xlabel('Distance')
    plt.ylabel('B Field Strength')
    plt.legend()
    plt.show()


def plot_sor_convergence():
    omega, times = read_from_file('sor_convergence.txt')

    min_point = np.argmin(times)
    minimum = round(omega[min_point], 3)

    plt.plot(omega, times)
    plt.xlabel('$\omega$')
    plt.ylabel('Time to convergence')
    plt.title('Time to convergence against $\omega$')
    plt.text(1.4, 150, 'The minimum $\omega$ is ' + str(minimum))
    plt.show()


plot_free_energy()
plot_poisson()
plot_magnetic()
plot_gs()
plot_sor_convergence()
