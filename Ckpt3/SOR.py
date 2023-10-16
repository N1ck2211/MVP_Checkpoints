import numpy as np
from ZweiFisch import Poisson


class SOR(Poisson):
    """
    Class to implement the SOR adaptation of Gauss-Seidel
    Extends Poisson
    """
    def __init__(self, size, nstep, omega, accuracy):
        Poisson.__init__(self, size, nstep, accuracy)
        self.omega = omega

    def nearest_neighbours_3D(self, grid, i, j, k):
        """
        Function to find the nearest neighbours of a given point in the array.

        :param grid: the array of points
        :param i: i coordinate
        :param j: j coordinate
        :param k: k coordinate
        :return: array of all 6 neighbours.
        """
        size = self.size

        up = grid[(i - 1) % size, j, k]
        down = grid[(i + 1) % size, j, k]
        left = grid[i, (j - 1) % size, k]
        right = grid[i, (j + 1) % size, k]
        forward = grid[i, j, (k + 1) % size]
        back = grid[i, j, (k - 1) % size]

        neighbours = np.array([up, down, left, right, forward, back])
        return neighbours

    def calc_potential(self):
        size = self.size
        omega = self.omega
        old_pot = np.sum(self.potential_grid)

        new_grid = self.potential_grid

        for i in range(size):
            for j in range(size):
                for k in range(size):
                    end = size - 1
                    if i == 0 or j == 0 or k == 0 or i == end or j == end or k == end:
                        self.potential_grid[i, j, k] = 0.0
                        self.charge_dist[i, j, k] = 0.0
                    else:
                        rho = self.charge_dist[i, j, k]

                        first_neighbours = self.nearest_neighbours_3D(self.potential_grid, i, j, k)
                        sum_one = first_neighbours[0] + first_neighbours[2] + first_neighbours[4]
                        second_neighbours = self.nearest_neighbours_3D(new_grid, i, j, k)
                        sum_two = second_neighbours[1] + second_neighbours[3] + second_neighbours[5]

                        neighbours_sum = sum_one + sum_two

                        phi_GS = 1 / 6 * (neighbours_sum + rho)

                        self.potential_grid[i, j, k] = omega * phi_GS + ((1 - omega) * self.potential_grid[i, j, k])

                        new_grid = self.potential_grid

        new_pot = np.sum(self.potential_grid)

        self.error = np.abs(new_pot - old_pot)

    def final_n(self):
        time = 0
        for n in range(self.nstep):
            self.calc_potential()

            print(self.error)

            if self.error < self.accuracy:
                print(n)
                time = n
                break

        return time


def convergence_test():
    omegas = np.arange(0.8, 2, 0.05)

    f = open('sor_convergence.txt', 'w')
    f.close()

    f = open('sor_convergence.txt', 'a')

    for omega in omegas:
        sor = SOR(20, 1000, omega, 0.01)
        sor.initialise()
        time = sor.final_n()
        f.write('%lf %lf \n' % (omega, time))

    f.close()

