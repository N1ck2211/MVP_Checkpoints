import matplotlib
import random
import numpy as np
import matplotlib.pyplot as plt

matplotlib.use('TKAgg')

J = 1.0


def graph_spins(size):
    spin = np.ones((size, size), dtype=float)
    return spin


def error_calc(sample_list):
    list_vals = sample_list
    squared_vals = sample_list ** 2

    mean = np.sum(list_vals) / len(list_vals)
    mean_squared = mean ** 2

    squared_vals_mean = np.sum(squared_vals) / len(squared_vals)

    return mean_squared, squared_vals_mean


class Glauber:
    def __init__(self, length, temp, nstep):
        self.len = length
        self.temp = temp
        self.nstep = nstep

    def init_spins(self):
        len_x = self.len
        len_y = self.len

        spin = np.zeros((len_x, len_y), dtype=float)

        # initialise spins randomly

        for i in range(len_x):
            for j in range(len_y):
                r = random.random()
                if r < 0.5:
                    spin[i, j] = -1
                if r >= 0.5:
                    spin[i, j] = 1

        return spin

    def find_energy(self, choice_spin, ispin, jspin, spin):
        len_x = self.len
        len_y = self.len

        # find the nearest neighbours
        down = spin[(ispin + 1) % len_y, jspin]
        up = spin[(ispin - 1) % len_y, jspin]
        right = spin[ispin, (jspin + 1) % len_x]
        left = spin[ispin, (jspin - 1) % len_x]

        # calculate ∆E
        energy = -J * ((up + down + left + right) * choice_spin)

        return energy

    def glauber(self, spin):
        # print("g")
        len_x = self.len
        len_y = self.len

        # choose a random spin
        itrial = np.random.randint(0, len_x)
        jtrial = np.random.randint(0, len_y)

        # take that spin and the flipped version
        init_spin = spin[itrial, jtrial]
        new_spin = -init_spin

        # calculate ∆E
        init_energy = self.find_energy(init_spin, itrial, jtrial, spin)
        new_energy = self.find_energy(new_spin, itrial, jtrial, spin)

        delt_E = (new_energy - init_energy)

        return delt_E, itrial, jtrial

    def metropolis(self, spin):
        delta_E, itrial, jtrial = self.glauber(spin)

        kT = self.temp

        rand = random.random()
        prob = np.exp((-delta_E) / kT)

        if delta_E <= 0:
            spin[itrial, jtrial] = spin[itrial, jtrial] * -1

        elif rand <= prob:
            spin[itrial, jtrial] = spin[itrial, jtrial] * -1

        else:
            spin[itrial, jtrial] = spin[itrial, jtrial]

        return spin

    def animate_glauber(self, spin):
        len_x = self.len
        len_y = self.len
        nstep = self.nstep

        for n in range(nstep):
            if n % 10 == 0:
                #       update measurements
                #       dump output
                f = open('spins.dat', 'w')
                for i in range(len_x):
                    for j in range(len_y):
                        spin = self.metropolis(spin)
                        f.write('%d %d %lf\n' % (i, j, spin[i, j]))
                f.close()
                #       show animation
                plt.cla()
                im = plt.imshow(spin, cmap='magma', animated=True)

                plt.draw()

                plt.pause(0.0001)

    def magnetisation(self, spin):
        # Function to calculate the magnetisation at a given temperature:
        len_x = self.len
        len_y = self.len
        nstep = self.nstep

        m_range = []

        # iterate through the first 100 sweeps without taking data
        for n in range(110):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

        # every 10 sweeps take a value for M and add to array
        for n in range(nstep):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

            if n % 10 == 0:
                m = np.sum(spin)
                m_range.append(np.abs(m))

        # calculate the average value of M, <M>:
        m_range = np.array(m_range)
        averageM = np.sum(m_range) / len(m_range)

        return averageM

    def susceptibility(self, spin):
        # Function to susceptibility the magnetisation at a given temperature:

        len_x = self.len
        len_y = self.len
        kT = self.temp
        nstep = self.nstep

        m_range = []
        squ_range = []

        # total spins, N:
        N = len_x * len_y

        # find value of initial fraction in suscep. equation
        term_one = 1 / (N*kT)

        # iterate through the first 100 steps without taking data
        for n in range(110):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

        # every 10 steps take a value for M and add to array
        for n in range(nstep):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

            if n % 10 == 0:
                M = np.sum(spin)
                M_sq = M ** 2

                m_range.append(np.abs(M))
                squ_range.append(M_sq)

        # calculate the average value of <M>, and <M^2>:
        m_range = np.asarray(m_range)
        squ_range = np.asarray(squ_range)

        average_m = np.sum(m_range) / len(m_range)
        average_sq = np.sum(squ_range) / len(squ_range)

        average_m_sq = average_m ** 2

        # Compute the second term in the equation then multiply w/ term one to yield suscep.
        term_two = average_sq - average_m_sq

        chi = term_one * term_two

        return chi

    def average_energy(self, spin):
        # Function to calculate the energy at a given temperature:
        len_x = self.len
        len_y = self.len
        nstep = self.nstep

        energy_range = []
        iter_energies = []

        # iterate through the first 100 steps without taking data
        for n in range(110):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

        # every 10 steps compute a value for energy by finding the sum of each spin's energy
        for n in range(nstep):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

            if n % 10 == 0:
                for i in range(len_x):
                    for j in range(len_y):
                        en_spin = self.find_energy(spin[i, j], i, j, spin)
                        iter_energies.append(en_spin)
                iter_energies = np.array(iter_energies)
                energy = np.sum(iter_energies)/2

                energy_range.append(energy)

            iter_energies = []

        # find the average energy at this temperature
        energy_range = np.asarray(energy_range)
        average_energy = np.sum(energy_range) / len(energy_range)

        return average_energy

    def heat_capacity(self, spin):
        # Function to calculate the heat capacity at a given temperature:
        len_x = self.len
        len_y = self.len

        kT = self.temp
        nstep = self.nstep

        energy_range = []
        iter_energies = []
        sq_energy_range = []

        term_one = 1 / (kT * kT)

        # iterate through the first 100 steps without taking data
        for n in range(110):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

        # every 10 steps compute a value for energy by finding the sum of each spin's energy

        for n in range(nstep):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.metropolis(spin)

            if n % 10 == 0:
                for i in range(len_x):
                    for j in range(len_y):
                        en_spin = self.find_energy(spin[i, j], i, j, spin)
                        iter_energies.append(en_spin)

                iter_energies = np.asarray(iter_energies)
                energy = np.sum(iter_energies) / 2
                energy_squared = energy ** 2
                energy_range.append(energy)
                sq_energy_range.append(energy_squared)

            iter_energies = []

        # calculate the average value of <E>, and <E^2>:
        energy_range = np.asarray(energy_range)
        sq_energy_range = np.asarray(sq_energy_range)

        average_energy = np.sum(energy_range) / len(energy_range)
        average_squ_energy = np.sum(sq_energy_range) / len(sq_energy_range)

        energy_sq = average_energy ** 2

        # Compute the second term in the equation then multiply w/ term one to yield heat capacity.
        term_two = average_squ_energy - energy_sq

        heat_c = term_one * term_two

        # Compute the error:
        sample_list = []
        c_list = []
        num_k = 1000

        for k in range(num_k):
            random_vals = np.random.randint(0, len(energy_range), size=len(energy_range))
            for i in random_vals:
                sample_list.append(energy_range[i])

            sample_list = np.array(sample_list)

            mean_squared, squared_vals_mean = error_calc(sample_list)

            c = term_one * (squared_vals_mean - mean_squared)
            c_list.append(c)

            sample_list = []

        c_list = np.array(c_list)
        c_sq_list = c_list ** 2

        av_c = np.sum(c_list) / len(c_list)
        av_c_sq_list = np.sum(c_sq_list) / len(c_sq_list)

        av_c_sq = av_c ** 2

        error = np.sqrt(av_c_sq_list - av_c_sq)

        return heat_c, error
