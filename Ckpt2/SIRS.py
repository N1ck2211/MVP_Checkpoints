import matplotlib
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap

matplotlib.use('TKAgg')

'''
Note: 0 means susceptible, 1 means infected, 2 means recovered 
'''

cmap = cm.get_cmap('gist_earth', 256)
new_colours = cmap(np.linspace(0, 1, 256))
vacc = np.array([(39 / 256, 236 / 256, 245 / 256, 0.8)])
new_colours[:-1, :] = vacc
newcmp = ListedColormap(new_colours)


class SIRS:
    def __init__(self, size, nstep, p1, p2, p3, im_frac):
        self.size = size
        self.nstep = nstep
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.im_frac = im_frac
        self.immune_coords = 0

    def init_cells_random(self):
        """Function to give a random set up of S, I, R"""

        size = self.size
        cell = np.zeros((size, size), dtype=float)

        # initialise cells randomly

        for i in range(size):
            for j in range(size):
                r = random.random()
                if r < 0.25:
                    cell[i, j] = 1
                if 0.25 < r < 0.3:
                    cell[i, j] = 2
                if r >= 0.3:
                    cell[i, j] = 0

        # one ill boi in a random place
        # i = np.random.randint(0, size)
        # j = np.random.randint(0, size)
        # cell[i, j] = 1

        return cell

    def init_cells_immune_frac(self):
        """
        Function to produce a randomly distributed set of SIR cells, but with a predetermined fraction remaining
        'immune'.
        :return: cell array
        """

        size = self.size
        N = self.size * self.size
        im_frac = self.im_frac

        cell = np.zeros((size, size), dtype=float)
        num_immune = int(N * im_frac)
        immune_list = []

        # generate the coordinates of immune cells:
        for im in range(0, num_immune):
            rand_i = random.randint(0, size - 1)
            rand_j = random.randint(0, size - 1)

            if [rand_i, rand_j] in immune_list:
                im -= 1
                continue
            else:
                immune_list.append([rand_i, rand_j])

        # create random array of SIR:
        for i in range(size):
            for j in range(size):
                r = random.random()
                if r < 0.3333:
                    cell[i, j] = 1
                if 0.3333 < r < 0.6666:
                    cell[i, j] = 2
                if r >= 0.6666:
                    cell[i, j] = 0

        # MAke every coordinate in the immune list an 'immune' cell
        for coord in immune_list:
            cell[coord[0], coord[1]] = 3

        self.immune_coords = immune_list

        return cell

    @staticmethod
    def nearest_neighbours(cell, i, j):
        """
        Function to return the values of a given cells nearest neighbours
        :param cell: cell grid
        :param i: i coordinate of target
        :param j: j coordinate of target
        :return: list of neighbours' values.
        """

        length = len(cell)

        # find 4 adjacent cells
        up = cell[(i - 1) % length, j]
        down = cell[(i + 1) % length, j]
        left = cell[i, (j - 1) % length]
        right = cell[i, (j + 1) % length]

        neighbours = np.array([up, down, left, right])
        return neighbours

    def sirs(self, cell):
        """
        Function to apply SIRS rules
        :param cell: cells grid
        :return: the cells grid after implementing SIRS rules
        """

        size = self.size

        # random place
        i = np.random.randint(0, size)
        j = np.random.randint(0, size)

        target = cell[i, j]

        neighbours = self.nearest_neighbours(cell, i, j)

        # print(target, neighbours)
        if target == 0:
            self.susceptible_rules(cell, i, j, neighbours)
        elif target == 1:
            self.infected_rules(cell, i, j)
        elif target == 2:
            self.recovered_rules(cell, i, j)

        return cell

    def sirs_with_immunity(self, cell):
        """
        Function to apply SIRS rules, modified to not change immune cells
        :param cell: cells grid
        :return: the cells grid after implementing SIRS rules
        """

        size = self.size
        im_coords = self.immune_coords

        # random place
        i = np.random.randint(0, size)
        j = np.random.randint(0, size)

        if [i, j] not in im_coords:

            target = cell[i, j]

            neighbours = self.nearest_neighbours(cell, i, j)

            # print(target, neighbours)
            if target == 0:
                self.susceptible_rules(cell, i, j, neighbours)
            elif target == 1:
                self.infected_rules(cell, i, j)
            elif target == 2:
                self.recovered_rules(cell, i, j)

        return cell

    def recovered_rules(self, cell, i, j):
        p3 = self.p3
        r3 = random.random()
        if r3 <= p3:
            cell[i, j] = 0

    def infected_rules(self, cell, i, j):
        p2 = self.p2
        r2 = random.random()
        if r2 <= p2:
            cell[i, j] = 2

    def susceptible_rules(self, cell, i, j, neighbours):
        p1 = self.p1
        r1 = random.random()
        if np.any(neighbours == 1) and r1 <= p1:
            cell[i, j] = 1

    def infected_fraction(self, cell):
        size = self.size
        N = size * size

        infected_list = []

        # iterate through the first 100 sweeps without taking data
        for n in range(100):
            for i in range(size):
                for j in range(size):
                    cell = self.sirs(cell)

        test_inf = np.sum(cell, where=(cell == 1))

        if test_inf == 0:
            average_infected_per_N = 0

        else:
            # take a value for I and add to array
            for n in range(self.nstep):
                for i in range(size):
                    for j in range(size):
                        cell = self.sirs(cell)

                infected = np.sum(cell, where=(cell == 1))
                infected_list.append(infected)

            infected_list = np.array(infected_list)
            average_infected = np.sum(infected_list) / len(infected_list)

            average_infected_per_N = average_infected / N

            print(average_infected_per_N)

        return average_infected_per_N

    def infected_fraction_with_immunity(self, cell):
        """
        Modifies the infected_fraction function to account for immune cells
        :param cell: cells gird
        :return: the average fraction of the grid that is infected
        """
        size = self.size
        N = size * size

        infected_list = []

        # iterate through the first 100 sweeps without taking data
        for n in range(100):
            for i in range(size):
                for j in range(size):
                    cell = self.sirs_with_immunity(cell)

        test_inf = np.sum(cell, where=(cell == 1))

        if test_inf == 0:
            average_infected_per_N = 0

        else:
            # take a value for I and add to array
            for n in range(self.nstep):
                for i in range(size):
                    for j in range(size):
                        cell = self.sirs_with_immunity(cell)

                infected = np.sum(cell, where=(cell == 1))
                infected_list.append(infected)

            infected_list = np.array(infected_list)

            average_infected = np.sum(infected_list) / len(infected_list)
            average_infected_per_N = average_infected / N

        return average_infected_per_N

    @staticmethod
    def error_calc(sample_list):
        """
        finds the <m>^2 and <m^2> of any list
        :param sample_list:
        :return: <m>^2 and <m^2>
        """
        list_vals = sample_list
        squared_vals = sample_list ** 2

        mean = np.sum(list_vals) / len(list_vals)
        mean_squared = mean ** 2

        squared_vals_mean = np.sum(squared_vals) / len(squared_vals)

        return mean_squared, squared_vals_mean

    def variance(self, cell):
        """
        Finds the variance and error in a cut of the phase
        :param cell: cells grid
        :return: variance, error
        """
        size = self.size

        # iterate through the first 100 sweeps without taking data
        for n in range(100):
            for i in range(size):
                for j in range(size):
                    cell = self.sirs(cell)

        test_inf = np.sum(cell, where=(cell == 1))

        # cuts calculation short if there are no infected after stabilisation
        if test_inf == 0:
            variance = 0
            error = 0

        else:
            error, variance = self.calculate_variance(cell)

        return variance, error

    def calculate_variance(self, cell):
        size = self.size
        N = size * size

        infected_list = []
        squared_infected_list = []

        # take a value for I and add to array
        for n in range(self.nstep):
            for i in range(size):
                for j in range(size):
                    cell = self.sirs(cell)

            # print(n)
            infected = np.sum(cell, where=(cell == 1))
            squared_infected = infected ** 2

            infected_list.append(infected)
            squared_infected_list.append(squared_infected)

        infected_list = np.array(infected_list)

        average_infected = np.sum(infected_list) / len(infected_list)
        average_infected_squared = average_infected ** 2

        squared_infected_list = np.array(squared_infected_list)
        average_squared_infected = np.sum(squared_infected_list) / len(squared_infected_list)

        variance = (average_squared_infected - average_infected_squared) / N

        error = self.variance_error(N, infected_list)

        return error, variance

    def variance_error(self, N, infected_list):
        # find the error with bootstrap:
        sample_list = []
        var_sample_list = []
        num_k = 1000

        for k in range(num_k):
            random_vals = np.random.randint(0, len(infected_list), size=len(infected_list))
            for i in random_vals:
                sample_list.append(infected_list[i])

            sample_list = np.array(sample_list)

            mean_squared, squared_vals_mean = self.error_calc(sample_list)

            var_sample = (squared_vals_mean - mean_squared) / N
            var_sample_list.append(var_sample)

            sample_list = []

        var_sample_list = np.array(var_sample_list)
        var_sq_list = var_sample_list ** 2

        av_var = (np.sum(var_sample_list) / len(var_sample_list))
        av_var_sq_list = np.sum(var_sq_list) / len(var_sq_list)

        av_var_sq = av_var ** 2
        error = np.sqrt(av_var_sq_list - av_var_sq)

        return error

    def animate(self, cell):
        """
        function to animate the standard SIRS cells
        :param cell: cells grid
        """

        size = self.size
        for n in range(self.nstep):

            for i in range(size):
                for j in range(size):
                    cell = self.sirs(cell)

            #       show animation
            plt.cla()
            im = plt.imshow(cell, cmap='gist_earth', animated=True)

            plt.draw()

            plt.pause(0.0001)

    def animate_immune(self, cell):
        """
        altered animate to show immune cells
        :param cell: cells grid
        """

        size = self.size
        for n in range(self.nstep):

            for i in range(size):
                for j in range(size):
                    cell = self.sirs_with_immunity(cell)

            #       show animation
            plt.cla()
            im = plt.imshow(cell, cmap='gist_earth', animated=True)
            plt.draw()

            plt.pause(0.0001)


def calc_phases():
    """
    iterate through a range of p1 and p3 values, calculating the infected fraction at each point
    """

    p1_range = np.arange(0, 1, 0.05)
    p3_range = np.arange(0, 1, 0.05)

    f = open('phases', 'w')
    f.close()

    f = open('phases', 'a')

    for p1 in p1_range:
        for p3 in p3_range:
            print(p1, p3)
            sirs = SIRS(50, p1, 0.5, p3, 0, 0.0)
            cells = sirs.init_cells_random()
            inf = sirs.infected_fraction(cells)
            f.write('%f %f %lf\n' % (p1, p3, inf))

    f.close()


def calc_variances():
    """
    iterate through a range of p1 values, calculating the variance and error at each point
    """
    p1s = np.arange(0.2, 0.51, 0.01)

    f = open('variances.txt', 'w')
    f.close()

    for p1 in p1s:
        print(p1)
        f = open('variances.txt', 'a')
        sirs = SIRS(50, p1, 0.5, 0.5, 0.5, 0)
        cells = sirs.init_cells_random()
        var, err = sirs.variance(cells)
        f.write('%f %f %f\n' % (p1, var, err))

        f.close()


def calc_infected_fraction_with_immune():
    """
    iterate through a range of immune cell fractions, calculating the infected fraction at each point
    """
    im_fractions = np.arange(0, 1.1, 0.01)

    f = open('infected_vs_immune_frac.txt', 'w')
    f.close()

    for frac in im_fractions:
        print(frac)
        f = open('infected_vs_immune_frac.txt', 'a')
        sirs = SIRS(50, 10000, 0.5, 0.5, 0.5, frac)
        cells = sirs.init_cells_immune_frac()
        inf_frac = sirs.infected_fraction_with_immunity(cells)
        f.write('%f %f\n' % (frac, inf_frac))

    f.close()


'''
Uncomment to calculate quantities and write to file.
'''
# calc_phases()
# calc_variances()
# calc_infected_fraction_with_immune()


'''
Absorbing state: 0.1, 0.5, 0.1
Dynamic Equilibrium: 1, 1, 1
Waves: 0.9, 0.6, 0.1
'''
