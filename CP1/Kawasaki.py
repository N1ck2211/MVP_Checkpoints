import matplotlib
import random
import numpy as np
import matplotlib.pyplot as plt

matplotlib.use('TKAgg')

J = 1.0


def kawasaki_graph_spins(size):
    # function to initialise the NxN array with half spin up, half down:

    # initialise array:
    spin = np.ones((size, size), dtype=float)

    # Find a random start value in the first half of the matrix:
    start_i = random.randint(0, size / 2)
    start_j = random.randint(0, size / 2)

    # calculate end values such that half the matrix is covered:
    end_i = (start_i + (size / 2)) % size
    end_j = (start_j + (size / 2)) % size

    state = 1

    # loop through the array, setting any spins between the start and end coordinates to -1
    for i in range(size):
        for j in range(size):
            if i == start_i and j == start_j:
                state = -1
            if i == end_i and j == end_j:
                state = 1

            spin[i, j] = state

    return spin


def error_calc(sample_list):
    # Function to find the <value>^2 and <value^2> of a given value
    list_vals = sample_list
    squared_vals = sample_list ** 2

    mean = np.sum(list_vals) / len(list_vals)
    mean_squared = mean ** 2

    squared_vals_mean = np.sum(squared_vals) / len(squared_vals)

    return mean_squared, squared_vals_mean


class Kawasaki:
    def __init__(self, length, temp, nstep):

        self.len = length
        self.temp = temp
        self.nstep = nstep

    def init_spins(self):
        #  initialise the spin array randomly; for animation
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

    def find_energy_spin(self, choice_spin, ispin, jspin, spins):
        # Find the energy at a given spin by summing over nearest neighbours
        len_x = self.len
        len_y = self.len

        down = spins[(ispin + 1) % len_y, jspin]
        up = spins[(ispin - 1) % len_y, jspin]
        right = spins[ispin, (jspin + 1) % len_x]
        left = spins[ispin, (jspin - 1) % len_x]

        # calculate E
        energy = -J * ((up + down + left + right) * choice_spin)

        return energy

    def get_neighbours_pos(self, ispin, jspin):
        # Finds the positions a spins nearest neighbours
        len_y = self.len
        len_x = self.len

        down = [(ispin + 1) % len_y, jspin]
        up = [(ispin - 1) % len_y, jspin]
        right = [ispin, (jspin + 1) % len_x]
        left = [ispin, (jspin - 1) % len_x]

        return up, down, left, right

    def check_for_degen(self, ispin_one, jspin_one, ispin_two, jspin_two):
        # Function checks if a spin is in the same position, or is a nearest neighbour as the partner spin
        alert = False

        up_one, down_one, left_one, right_one = self.get_neighbours_pos(ispin_one, jspin_one)

        if [ispin_one, jspin_one] == [ispin_two, jspin_two]:
            alert = True

        elif up_one == [ispin_two, jspin_two]:
            alert = True

        elif down_one == [ispin_two, jspin_two]:
            alert = True

        elif left_one == [ispin_two, jspin_two]:
            alert = True

        elif right_one == [ispin_two, jspin_two]:
            alert = True

        return alert

    def find_delta_e(self, spins):
        # Function to find the change in energy when two spins are flipped:
        len_x = self.len
        len_y = self.len
        correction = 0

        # choose a first random spin
        ispin_one = np.random.randint(0, len_x)
        jspin_one = np.random.randint(0, len_y)
        spin_one = spins[ispin_one, jspin_one]

        # choose a second random spin
        ispin_two = np.random.randint(0, len_x)
        jspin_two = np.random.randint(0, len_y)
        spin_two = spins[ispin_two, jspin_two]

        while spin_one == spin_two:
            # choose a different spin:
            # choose a first random spin
            ispin_one = np.random.randint(0, len_x)
            jspin_one = np.random.randint(0, len_y)
            spin_one = spins[ispin_one, jspin_one]

            # choose a second random spin
            ispin_two = np.random.randint(0, len_x)
            jspin_two = np.random.randint(0, len_y)
            spin_two = spins[ispin_two, jspin_two]

        # check if neighbours or identical:
        alert = self.check_for_degen(ispin_one, jspin_one, ispin_two, jspin_two)

        if alert:
            correction = -4

        # find the initial energies of the spins:
        energy_spin_one_init = self.find_energy_spin(spin_one, ispin_one, jspin_one, spins)
        energy_spin_two_init = self.find_energy_spin(spin_two, ispin_two, jspin_two, spins)

        # find the energies with the spins flipped:
        energy_spin_one_new = self.find_energy_spin(spin_two, ispin_one, jspin_one, spins)
        energy_spin_two_new = self.find_energy_spin(spin_one, ispin_two, jspin_two, spins)

        # find deltaE for both spins and sum them
        delta_E_one = energy_spin_one_new - energy_spin_one_init
        delta_E_two = energy_spin_two_new - energy_spin_two_init

        delta_E = delta_E_one + delta_E_two - correction

        return delta_E, ispin_one, jspin_one, ispin_two, jspin_two

    def kawasaki(self, spins):
        # Function to determine if a spin swap is carried out:
        delta_E, i_one, j_one, i_two, j_two = self.find_delta_e(spins)

        # set the initial states of the spins
        spin_one = spins[i_one, j_one]
        spin_two = spins[i_two, j_two]

        kT = self.temp

        # set the random number and the probability it will be compared to
        rand = random.random()
        prob = np.exp((-delta_E) / kT)

        # if the conditions are met, flip the spins, else keep them in their original place
        if delta_E <= 0:
            spins[i_one, j_one] = spin_two
            spins[i_two, j_two] = spin_one

        elif rand <= prob:
            spins[i_one, j_one] = spin_two
            spins[i_two, j_two] = spin_one

        else:
            spins[i_one, j_one] = spin_one
            spins[i_two, j_two] = spin_two

        return spins

    def animate(self, spin):
        # function to display an animation of the ising model:

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
                        spin = self.kawasaki(spin)
                        f.write('%d %d %lf\n' % (i, j, spin[i, j]))
                f.close()
                #       show animation
                plt.cla()
                im = plt.imshow(spin, cmap='magma', animated=True)
                plt.draw()
                plt.pause(0.0001)

    def kawasaki_mag(self, spin):
        # Function to calculate the magnetisation at a given temperature:
        len_x = self.len
        len_y = self.len
        nstep = self.nstep

        M_range = []
        # iterate through the first 100 sweeps without taking data
        for n in range(110):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.kawasaki(spin)

        # every 10 sweeps take a value for M and add to array
        for n in range(nstep):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.kawasaki(spin)

            if n % 10 == 0:
                M = np.sum(spin)
                M_range.append(np.abs(M))

        # calculate the average value of M, <M>:
        M_range = np.array(M_range)
        averageM = np.sum(M_range) / len(M_range)

        return averageM

    def average_energy_kawasaki(self, spin):
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
                    spin = self.kawasaki(spin)

        # every 10 steps compute a value for energy by finding the sum of each spin's energy
        for n in range(nstep):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.kawasaki(spin)

            if n % 10 == 0:
                for i in range(len_x):
                    for j in range(len_y):
                        en_spin = self.find_energy_spin(spin[i, j], i, j, spin)
                        iter_energies.append(en_spin)
                iter_energies = np.array(iter_energies)

                # add the computed energy to energy range
                energy = np.sum(iter_energies) / 2
                energy_range.append(energy)

            iter_energies = []

        # find the average energy at this temperature
        energy_range = np.asarray(energy_range)
        average_energy = np.sum(energy_range) / len(energy_range)

        return average_energy

    def heat_capacity_kawasaki(self, spin):
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
                    spin = self.kawasaki(spin)

        # every 10 steps compute a value for energy by finding the sum of each spin's energy

        for n in range(nstep):
            for i in range(len_x):
                for j in range(len_y):
                    spin = self.kawasaki(spin)

            if n % 10 == 0:
                for i in range(len_x):
                    for j in range(len_y):
                        en_spin = self.find_energy_spin(spin[i, j], i, j, spin)
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

