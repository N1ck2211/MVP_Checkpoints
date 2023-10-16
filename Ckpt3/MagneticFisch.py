import numpy as np
from ZweiFisch import Poisson
import random
import matplotlib.pyplot as plt


class MagneticPoisson(Poisson):
    """
    Class to use the poisson equation to find the potential and magnetic field for a wire through the z axis.
    Extends Poisson.
    """
    def __init__(self, size, nstep, accuracy):
        Poisson.__init__(self, size, nstep, accuracy)

    def initialise(self):
        """
        Function to initialise the system with a wire through the z axis.
        """
        size = self.size
        midpoint = int(self.size / 2)

        for i in range(size):
            for j in range(size):
                for k in range(size):
                    end = size - 1
                    if j == 0 or k == 0 or j == end or k == end:
                        self.potential_grid[i, j, k] = 0.0
                        self.charge_dist[i, j, k] = 0.0

                    else:
                        noise = random.uniform(-0.01, 0.01)
                        self.potential_grid[i, j, k] = noise

        self.charge_dist[:, midpoint, midpoint] = 1.0

    def endpoint(self):
        """
         Determine when the threshold for convergence has been reached and display plots/read to file.
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

        plt.contourf(x, y, self.potential_grid[int(size / 2)], 40, cmap='magma')
        plt.colorbar()
        plt.show()

        Bx, By = np.gradient(self.potential_grid[1, :, :])

        B = np.sqrt((Bx**2) + (By**2))

        X, Y = np.meshgrid(x, y)
        u = -Bx / B
        v = By / B

        plt.quiver(X, Y, u, v)
        plt.show()

        name = 'B'
        self.distance_strength(midpoint, B, name)



