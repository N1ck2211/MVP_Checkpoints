import numpy as np
import matplotlib.pyplot as plt
from ZweiFisch import Poisson


class GaussSeidel(Poisson):
    """
    Class to calculate the potential and electric field using the Gauss-Seidel algorithm.
    Extends Poisson.
    """

    def __init__(self, size, nstep, accuracy):
        Poisson.__init__(self, size, nstep, accuracy)

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
        """
        Apply the Gauss-Seidel algorithm to calculate the potential
        :return:
        """
        size = self.size
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
                        self.potential_grid[i, j, k] = 1 / 6 * (neighbours_sum + rho)

                        new_grid = self.potential_grid

        new_pot = np.sum(self.potential_grid)

        # Find the difference between this step and the last; needed for future functions
        self.error = np.abs(new_pot-old_pot)

    def endpoint(self):
        size = self.size

        x = np.arange(0, size, 1)
        y = np.arange(0, size, 1)

        midpoint = int(size / 2)

        for n in range(self.nstep):
            self.calc_potential()

            print(self.error)

            if self.error < self.accuracy:
                print(n)
                break

        plt.contourf(x, y, self.potential_grid[midpoint], 40, cmap='magma')
        plt.colorbar()
        plt.show()

        Ex = np.gradient(self.potential_grid, axis=2)
        Ey = np.gradient(self.potential_grid, axis=1)

        Emag = np.sqrt(Ex[int(midpoint)] ** 2 + Ey[int(midpoint)] ** 2)

        X, Y = np.meshgrid(x, y)

        u = -Ex[int(midpoint)] / Emag
        v = -Ey[int(midpoint)] / Emag

        plt.quiver(X, Y, u, v)
        plt.show()

        E = np.sqrt((Ex[midpoint] ** 2) + (Ey[midpoint] ** 2))
        name = 'G'

        self.distance_strength(midpoint, size, E, name)
