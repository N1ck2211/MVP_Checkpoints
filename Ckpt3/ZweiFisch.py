import numpy as np
import random
import matplotlib.pyplot as plt


class Poisson:
    """
    Class to implement the Poisson algorith to calculate the potential and electric field from a charge distribution.
    forms the basis from which the other boundary condition implementations are extended from.
    """

    def __init__(self, size, nstep, accuracy):
        self.error = None
        self.size = size
        self.nstep = nstep
        self.dt = 1 / 2
        self.dx = 1
        self.epsilon = 1
        self.potential_grid = np.zeros((size, size, size), dtype=float)
        self.charge_dist = np.zeros((size, size, size), dtype=float)
        self.e_field = np.zeros((size, size, size), dtype=float)
        self.accuracy = accuracy

    # Using Dirichlet?
    def initialise(self):
        """
        Function to set up the 3d array with either a single charge in the centre, or with a random set of charges
        """
        size = self.size

        # place a charge in the centre
        midpoint = int(self.size / 2)
        self.charge_dist[midpoint, midpoint, midpoint] = 1

        for i in range(size):
            for j in range(size):
                for k in range(size):
                    end = size - 1

                    # implements the boundary conditions:
                    if i == 0 or j == 0 or k == 0 or i == end or j == end or k == end:
                        self.potential_grid[i, j, k] = 0.0
                        self.charge_dist[i, j, k] = 0.0

                    else:
                        noise = random.uniform(-0.01, 0.01)
                        self.potential_grid[i, j, k] = noise
                        # random_i = random.randint(1, size - 1)
                        # random_j = random.randint(1, size - 1)
                        # random_k = random.randint(1, size - 1)
                        # random_place = random.random()
                        # #
                        # if random_place <= 0.02:
                        #     self.charge_dist[random_i, random_j, random_k] = 1

    def calc_potential(self):
        """
        uses the roll function to implement the poisson equation
        """

        size = self.size
        old_pot = np.sum(self.potential_grid)

        # Find the neighbours of each element in the grid:
        forward = np.roll(self.potential_grid, 1, axis=1)
        back = np.roll(self.potential_grid, -1, axis=1)
        left = np.roll(self.potential_grid, 1, axis=2)
        right = np.roll(self.potential_grid, -1, axis=2)
        up = np.roll(self.potential_grid, 1, axis=0)
        down = np.roll(self.potential_grid, -1, axis=0)

        neighbour_sum = up + down + left + right + forward + back

        # Implement poisson
        self.potential_grid = 1 / 6 * (neighbour_sum + self.charge_dist)

        # Apply boundary conditions.
        self.potential_grid[::size - 1] = 0
        self.potential_grid[0::, ::size - 1] = 0
        self.potential_grid[0::, 0::, ::size - 1] = 0

        new_pot = np.sum(self.potential_grid)

        # Find the difference between this step and the last; needed for future functions
        self.error = np.abs(new_pot-old_pot)

    def animate(self):
        """
        Animate the behaviour of a slice through the mid-plane:
        """
        size = self.size
        for n in range(self.nstep):
            if n % 100 == 0:
                self.calc_potential()

                plt.cla()
                im = plt.imshow(self.potential_grid[int(size / 2)], cmap='magma', animated=True)
                plt.draw()
                plt.pause(0.0001)

    def endpoint(self):
        """
        Function determines when the error has reached the specified threshold. Call functions to show a contour plot
        of the potential, a vector plot of the field, and the data for potential and field at a distance from the centre
        """
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

        self.show_contour(midpoint, x, y)

        E, name = self.display_field(midpoint, x, y)

        self.distance_strength(midpoint, E, name)

    def display_field(self, midpoint, x, y):
        """
        Function to calculate and display the electric field as a vector plot

        :param midpoint: The centre of the array
        :param x: An array of the x-axis values
        :param y: An array of the y-axis values
        :return: E, an array of the E-field magnitudes, the name to tell the plotter which field to plot
        """

        # Use the gradient function to find the Electric field from the potential.
        Ex = np.gradient(self.potential_grid, axis=2)
        Ey = np.gradient(self.potential_grid, axis=1)

        E = np.sqrt((Ex[midpoint] ** 2) + (Ey[midpoint] ** 2))

        X, Y = np.meshgrid(x, y)

        u = -Ex[int(midpoint)] / E
        v = -Ey[int(midpoint)] / E

        plt.quiver(X, Y, u, v)
        plt.show()

        name = 'E'

        return E, name

    def show_contour(self, midpoint, x, y):
        plt.contourf(x, y, self.potential_grid[midpoint], 40, cmap='magma')
        plt.colorbar()
        plt.show()

    def distance_strength(self, midpoint, field, name):
        """
        Function to calculate the distance of a point to the midpoint, then write this and the field/potential to that
         point to a file.

        :param midpoint: the centre of the array, where the charge is located.
        :param field: the potential/field to be written to file
        :param name: the field/algorithm being used.
        """

        counter = 0

        # Writing the poisson potential and electric field data to file
        if name == 'E':
            pp = open('poisson_potential.txt', 'w')
            pp.close()

            f = open('electric_field_poisson.txt', 'w')

            f.close()

            pp = open('poisson_potential.txt', 'a')
            f = open('electric_field_poisson.txt', 'a')

            for i in range(len(self.potential_grid[midpoint])):
                for j in range(len(self.potential_grid[midpoint])):
                    px = midpoint - i
                    py = midpoint - j
                    r = np.sqrt((px ** 2) + (py ** 2))

                    pot_magnitude = self.potential_grid[midpoint][i][j]
                    field_magnitude = field[i][j]

                    pp.write('%lf %lf\n' % (r, pot_magnitude))
                    f.write('%lf %lf\n' % (r, field_magnitude))

                    counter += 1

            pp.close()
            f.close()

        # Writing the Gauss-Seidel potential and electric field data to file:
        elif name == 'G':
            pp = open('gs_potential.txt', 'w')
            pp.close()

            f = open('electric_field_gs.txt', 'w')
            f.close()

            pp = open('gs_potential.txt', 'a')
            f = open('electric_field_gs.txt', 'a')

            for i in range(len(self.potential_grid[midpoint])):
                for j in range(len(self.potential_grid[midpoint])):
                    px = midpoint - i
                    py = midpoint - j
                    r = np.sqrt((px ** 2) + (py ** 2))

                    pot_magnitude = self.potential_grid[midpoint][i][j]
                    field_magnitude = field[i][j]

                    pp.write('%lf %lf\n' % (r, pot_magnitude))
                    f.write('%lf %lf\n' % (r, field_magnitude))

                    counter += 1

            pp.close()
            f.close()

        # Writing the Poisson potential and magnetic field for a wire to file
        elif name == 'B':
            pp = open('magnetic_potential.txt', 'w')
            pp.close()

            f = open('magnetic_field.txt', 'w')
            f.close()

            pp = open('magnetic_potential.txt', 'a')
            f = open('magnetic_field.txt', 'a')

            for i in range(len(self.potential_grid[midpoint])):
                for j in range(len(self.potential_grid[midpoint])):
                    px = midpoint - i
                    py = midpoint - j
                    r = np.sqrt((px ** 2) + (py ** 2))

                    pot_magnitude = self.potential_grid[midpoint][i][j]
                    field_magnitude = field[i][j]

                    pp.write('%lf %lf\n' % (r, pot_magnitude))
                    f.write('%lf %lf\n' % (r, field_magnitude))

                    counter += 1

            pp.close()
            f.close()

        else:
            print('Something went wrong with the field name')

