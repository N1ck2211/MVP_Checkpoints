import matplotlib.pyplot as plt
import numpy as np
from progress.bar import Bar

from Kawasaki import Kawasaki
from Kawasaki import kawasaki_graph_spins
from Glauber import Glauber
from Glauber import graph_spins


def display():
    # Function to take in input and display the animation based on these parameters
    size = int(input("Input size n of nxn grid: "))
    temp = float(input("Input temperature: "))
    nstep = int(input('How many n steps?'))

    glaub = Glauber(size, temp, nstep)
    kawasaki = Kawasaki(size, temp, nstep)
    spins = glaub.init_spins()

    algo_choice = input('Using Glauber or Kawasaki dynamics? (G/K) ')

    if algo_choice == 'g' or algo_choice == 'G':
        glaub.animate_glauber(spins)
    elif algo_choice == 'k' or algo_choice == 'K':
        kawasaki.animate(spins)
    else:
        print("Invalid input, please input G or K ")


def graph_glauber():
    # Produces a graph and writes to the datafile based on the selected quantity
    size = int(input("Input size n of nxn grid: "))
    step = float(input("Input temperature step size: "))
    quantity = input("Which quantity would you like to graph? (M/S/E/C) ")
    nstep = int(input('How many n steps? '))

    t_range = np.arange(1, 3, step)

    m_vals = []
    chi_vals = []
    energy_vals = []

    c_vals = []
    c_errors = []

    if quantity == 'M' or quantity == 'm':
        bar = Bar('running magnetisation', fill='█', max=len(t_range))
        spins = graph_spins(size)

        f = open('magnetisation_glauber.txt', 'w')
        f.write('%s %s\n' % ('Temperature', 'Magnetisation'))
        f.close()

        for temp in t_range:
            f = open('magnetisation_glauber.txt', 'a')
            glaub = Glauber(size, temp, nstep)
            m_val = glaub.magnetisation(spins)
            m_vals.append(m_val)
            f.write('%f %f\n' % (temp, m_val))
            f.close()
            bar.next()
        bar.finish()

        plt.plot(t_range, m_vals, marker='o', label='Magnetisation')
        plt.title('Magnetisation against Temperature, \n Calculated from the Ising Model using Glauber Dynamics')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Magnetisation")

        plt.legend()
        plt.show()

    elif quantity == 'S' or quantity == 's':
        bar = Bar('running susceptibility', fill='█', max=len(t_range))
        spins = graph_spins(size)

        f = open('susceptibility_glauber.txt', 'w')
        f.write('%s %s\n' % ('Temperature', 'Susceptibility'))
        f.close()

        for temp in t_range:
            f = open('susceptibility_glauber.txt', 'a')
            glaub = Glauber(size, temp, nstep)
            chi_val = glaub.susceptibility(spins)
            chi_vals.append(chi_val)
            f.write('%f %f\n' % (temp, chi_val))
            f.close()
            bar.next()
        bar.finish()

        plt.plot(t_range, chi_vals, label='Susceptibility')
        plt.title('Susceptibility against Temperature, \n Calculated from the Ising Model using Glauber Dynamics')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Susceptibility")

        plt.legend()
        plt.show()

    elif quantity == 'E' or quantity == 'e':
        bar = Bar('running energy ', fill='█', max=len(t_range))
        spins = graph_spins(size)

        f = open('energy_glauber.txt', 'w')
        f.write('%s %s\n' % ('Temperature', 'Energy'))
        f.close()

        for temp in t_range:
            f = open('energy_glauber.txt', 'a')
            glaub = Glauber(size, temp, nstep)
            energy_val = glaub.average_energy(spins)
            energy_vals.append(energy_val)
            f.write('%f %f\n' % (temp, energy_val))
            f.close()
            bar.next()
        bar.finish()

        plt.plot(t_range, energy_vals, label='Energy')
        plt.title('Energy against Temperature, \n Calculated from the Ising Model using Glauber Dynamics')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Energy")

        plt.legend()
        plt.show()

    elif quantity == 'C' or quantity == 'c':
        bar = Bar('running heat capacity', fill='█', max=len(t_range))
        spins = graph_spins(size)

        f = open('capacity_glauber.txt', 'w')
        f.write('%s %s %s\n' % ('Temperature', 'Heat Capacity', 'Error'))
        f.close()

        for temp in t_range:
            f = open('capacity_glauber.txt', 'a')
            glaub = Glauber(size, temp, nstep)
            c_val, c_err = glaub.heat_capacity(spins)
            c_vals.append(c_val)
            c_errors.append(c_err)
            f.write('%f %f %f\n' % (temp, c_val, c_err))
            f.close()
            bar.next()
        bar.finish()

        plt.errorbar(t_range, c_vals, yerr=c_errors, ecolor='black', capsize=5, label='Heat Capacity')
        plt.title('Heat Capacity against Temperature, \n Calculated from the Ising Model using Glauber Dynamics')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Heat Capacity")

        plt.legend()
        plt.show()

    else:
        print("Not a valid quantity input")


def graph_kawasaki():
    size = int(input("Input size n of nxn grid: "))
    step = float(input("Input temperature step size: "))
    kawasaki_quantity = input("Which quantity would you like to graph? (M/E/C) ")
    nstep = int(input('How many n steps? '))

    t_range = np.arange(1, 3, step)
    m_vals = []
    energy_vals = []
    c_vals = []
    c_errors = []

    if kawasaki_quantity == 'M' or kawasaki_quantity == 'm':
        bar = Bar('running magnetisation', fill='█', max=len(t_range))
        spins = kawasaki_graph_spins(size)

        f = open('magnetisation_kawasaki.txt', 'w')
        f.write('%s %s\n' % ('Temperature', 'Magnetisation'))
        f.close()

        for temp in t_range:
            f = open('magnetisation_kawasaki.txt', 'a')
            kawasaki = Kawasaki(size, temp, nstep)
            m_val = kawasaki.kawasaki_mag(spins)
            m_vals.append(m_val)
            f.write('%f %f\n' % (temp, m_val))
            f.close()
            bar.next()
        bar.finish()

        plt.plot(t_range, m_vals, label='Magnetisation')
        plt.title('Magnetisation against Temperature, \n Calculated from the Ising Model using Kawasaki Dynamics')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Magnetisation")

        plt.legend()
        plt.show()

    elif kawasaki_quantity == 'E' or kawasaki_quantity == 'e':
        bar = Bar('running energy', fill='█', max=len(t_range))
        spins = kawasaki_graph_spins(size)
        f = open('energy_kawasaki.txt', 'w')
        f.write('%s %s\n' % ('Temperature', 'Energy'))
        f.close()

        for temp in t_range:
            f = open('energy_kawasaki.txt', 'a')
            kawasaki = Kawasaki(size, temp, nstep)
            energy_val = kawasaki.average_energy_kawasaki(spins)
            energy_vals.append(energy_val)
            f.write('%f %f\n' % (temp, energy_val))
            f.close()
            bar.next()
        bar.finish()

        plt.plot(t_range, energy_vals, label='Energy')
        plt.title('Energy against Temperature, \n Calculated from the Ising Model using Kawasaki Dynamics')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Energy")

        plt.legend()
        plt.show()

    elif kawasaki_quantity == 'C' or kawasaki_quantity == 'c':
        bar = Bar('running heat capacity', fill='█', max=len(t_range))
        spins = kawasaki_graph_spins(size)

        f = open('capacity_kawasaki.txt', 'w')
        f.write('%s %s %s\n' % ('Temperature', 'Heat Capacity', 'Error'))
        f.close()

        for temp in t_range:
            f = open('capacity_kawasaki.txt', 'a')
            kawasaki = Kawasaki(size, temp, nstep)
            c_val, c_err = kawasaki.heat_capacity_kawasaki(spins)
            c_vals.append(c_val)
            c_errors.append(c_err)
            f.write('%f %f %f\n' % (temp, c_val, c_err))
            f.close()
            bar.next()
        bar.finish()

        # plt.plot(t_range, c_vals, marker='o', label='Heat Capacity')
        plt.errorbar(t_range, c_vals, yerr=c_errors, ecolor='black', capsize=5, label='Heat Capacity')
        plt.title('Heat Capacity against Temperature, \n Calculated from the Ising Model using Kawasaki Dynamics')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Heat Capacity")

        plt.legend()
        plt.show()

    else:
        print("Not a valid quantity input")


def run():
    # Function to run the code - it takes in user input to determine which function to perform
    task = input('Do you want an animation of the ising model or a graphical output? (A/G) ')

    if task == 'A' or task == 'a':
        display()

    elif task == 'G' or task == 'g':
        dyn = input('Using Glauber or Kawasaki dynamics? (G/K) ')

        if dyn == 'G' or dyn == 'g':
            graph_glauber()
        elif dyn == 'K' or dyn == 'k':
            graph_kawasaki()
        else:
            print('Incorrect input, please input G or K')

    else:
        print('Incorrect input, please input A or G')


def glauber_graph_from_file(file):
    # function to produce graphs from a given datafile

    f = open(file)
    data = f.readlines()
    titles = data[0]
    data = data[1::]
    temp_list = []
    data_list = []
    error_list = []

    for i in range(20):
        temp_list.append(float(((data[i].split())[0])))
        data_list.append(float(((data[i].split())[1])))

    # checks to see if error-bars need to be added
    if len(data[0].split()) == 3:
        for i in range(20):
            error_list.append(float(((data[i].split())[2])))

    quantity = (titles.split())[1]

    if len(data[0].split()) == 3:
        plt.errorbar(temp_list, data_list, yerr=error_list, ecolor='black', capsize=5, label='Heat Capacity')
        plt.title('Heat Capacity against Temperature, \n Calculated from the Ising Model using Glauber Dynamics')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Heat Capacity')

    else:
        plt.plot(temp_list, data_list, marker='x', label=quantity)
        plt.title(quantity + ' against Temperature, \n Calculated from the Ising Model using Glauber Dynamics')
        plt.xlabel('Temperature (K)')
        plt.ylabel(quantity)

    plt.legend()
    plt.show()


def kawasaki_graph_from_file(file):
    # function to produce graphs from a given datafile

    f = open(file)
    data = f.readlines()
    titles = data[0]
    data = data[1::]
    temp_list = []
    data_list = []
    error_list = []

    for i in range(20):
        temp_list.append(float(((data[i].split())[0])))
        data_list.append(float(((data[i].split())[1])))

    # checks to see if error-bars need to be added
    if len(data[0].split()) == 3:
        for i in range(20):
            error_list.append(float(((data[i].split())[2])))

    quantity = (titles.split())[1]

    if len(data[0].split()) == 3:
        plt.errorbar(temp_list, data_list, yerr=error_list, ecolor='black', capsize=5, label='Heat Capacity')
        plt.title('Heat Capacity against Temperature, \n Calculated from the Ising Model using Kawasaki Dynamics')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Heat Capacity')

    else:
        plt.plot(temp_list, data_list, marker='x', label=quantity)
        plt.title(quantity + ' against Temperature, \n Calculated from the Ising Model using Kawasaki Dynamics')
        plt.xlabel('Temperature (K)')
        plt.ylabel(quantity)

    plt.legend()
    plt.show()


""" comment out run() if you want to plot graphs from the data file """

# run()

"""comment out graph from files functions if you want to run an animation or create new data"""

glauber_graph_from_file('magnetisation_glauber.txt')
glauber_graph_from_file('susceptibility_glauber.txt')
glauber_graph_from_file('energy_glauber.txt')
glauber_graph_from_file('capacity_glauber.txt')

kawasaki_graph_from_file('energy_kawasaki.txt')
kawasaki_graph_from_file('capacity_kawasaki.txt')




