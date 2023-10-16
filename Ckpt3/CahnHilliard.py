import numpy as np
import random
import matplotlib.pyplot as plt


class CahnHillard:
    """
    Class to implement the Cahn-Hilliard equation
    """
    def __init__(self, size, nstep, psi_o):
        self.size = size
        self.nstep = nstep
        self.psi_o = psi_o
        self.M = 0.1
        self.a = 0.1
        self.k = 0.1
        self.del_x = 1
        self.del_t = 1
        self.grid = np.zeros((size, size), dtype=float)
        self.mu_grid = np.zeros((size, size), dtype=float)

    def initialise_system(self):
        """
        Function to initialise the grid based on the inputted phi_o and some random noise
        """
        size = self.size

        for i in range(size):
            for j in range(size):
                noise = random.uniform(-0.01, 0.01)
                self.grid[i, j] = self.psi_o + noise

    def chemical_potential(self):
        """
        Function to calculate the array of chemical potential
        """
        a = self.a
        k = self.k
        del_x = self.del_x

        psi = self.grid

        # use roll to produce matrices of the neighbours
        left = np.roll(psi, (0, 1))
        right = np.roll(psi, (0, -1))
        up = np.roll(psi, 1, axis=0)
        down = np.roll(psi, -1, axis=0)

        neighbours_sum = up + down + left + right

        # calculate chemical potential
        self.mu_grid = (-a * psi + (a * (psi ** 3)) - (k / (del_x * del_x)) * (neighbours_sum - (4 * psi)))

    def update_psi(self):
        """
        Function to calculate the order parameter due to the chemical potential
        """
        M = self.M
        del_t = self.del_t
        del_x = self.del_x
        mu_grid = self.mu_grid

        # use roll to produce matrices of the neighbours
        left = np.roll(mu_grid, (0, 1))
        right = np.roll(mu_grid, (0, -1))
        up = np.roll(mu_grid, 1, axis=0)
        down = np.roll(mu_grid, -1, axis=0)

        mu_sum = up + down + left + right

        # Calculate the order parameter
        self.grid = self.grid + (M * del_t / (del_x * del_x)) * (mu_sum + (-4 * mu_grid))

    def free_energy(self):
        """
        Function to calculate the free energy
        """
        a = self.a
        k = self.k

        psi = self.grid

        # use roll to produce matrices of the neighbours
        left = np.roll(psi, (0, 1))
        right = np.roll(psi, (0, -1))
        up = np.roll(psi, 1, axis=0)
        down = np.roll(psi, -1, axis=0)

        neighbours_sum = up - down + left - right

        # calculate free energy
        f = -(a / 2 * (psi ** 2)) + (a / 4 * (psi ** 4)) + (k / 2 * (neighbours_sum / 2 * self.del_x) ** 2)
        f = np.sum(f)

        return f

    def animate(self):
        """
        Function to show an animation of the Cahn-Hilliard behaviour
        :return:
        """
        size = self.size

        for n in range(self.nstep):

            for i in range(0, size ** 2):
                self.chemical_potential()
                self.update_psi()

            # show animation
            plt.cla()
            im = plt.imshow(self.grid, cmap='BrBG', animated=True)
            plt.draw()

            plt.pause(0.0001)

    def get_f(self):
        """
        Function to find the free energy at a given step and write to file
        """
        nstep = self.nstep

        f_range = np.zeros(nstep)

        f = open('free_energy.txt', 'w')
        f.close()
        f = open('free_energy.txt', 'a')

        for n in range(nstep):
            print(n)
            self.chemical_potential()
            self.update_psi()
            free_energy = float(self.free_energy())
            f_range[n] = free_energy

            f.write('%lf %lf \n' % (n, free_energy))

        f.close()
