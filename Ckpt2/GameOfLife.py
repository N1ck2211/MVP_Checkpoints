import matplotlib
import numpy as np
import random
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')


class GameOfLife:
    def __init__(self, size, nstep):
        self.nstep = nstep
        self.size = size

    def init_cells_random(self):
        """
        :return: An array of randomly distributed alive and dead cells
        """
        size = self.size

        cell = np.zeros((size, size), dtype=float)

        # initialise cells randomly

        for i in range(size):
            for j in range(size):
                r = random.random()
                if r < 0.5:
                    cell[i, j] = -1
                if r >= 0.5:
                    cell[i, j] = 1

        return cell

    def init_cells_blinker(self):
        """
        :return: An array with only one blinker
        """
        size = self.size

        cell = np.ones((size, size), dtype=float)
        cell *= -1

        blink_pos = random.randint(5, size - 5)

        cell[blink_pos % size, blink_pos % size] = 1
        cell[blink_pos % size, (blink_pos + 1) % size] = 1
        cell[blink_pos % size, (blink_pos - 1) % size] = 1

        return cell

    def init_cells_glider(self):
        """
        :return: An array with only a single glider
        """
        size = self.size

        cell = np.ones((size, size), dtype=float)
        cell *= -1

        blink_pos = random.randint(1, size - 1)

        cell[blink_pos % size, blink_pos % size] = 1
        cell[blink_pos % size, (blink_pos + 1) % size] = 1
        cell[blink_pos % size, (blink_pos - 1) % size] = 1
        cell[(blink_pos - 2) % size, blink_pos % size] = 1
        cell[(blink_pos - 1) % size, (blink_pos + 1) % size] = 1

        return cell

    @staticmethod
    def neighbours(cell, i, j):
        """
        Function to find the value at all 8 of a cells neighbours
        :param cell: the cell grid
        :param i: x position
        :param j: y position
        :return: array of neighbour values
        """
        length = len(cell)

        # find all 8 neighbors
        up = cell[(i - 1) % length, j]
        down = cell[(i + 1) % length, j]
        left = cell[i, (j - 1) % length]
        right = cell[i, (j + 1) % length]

        l_up = cell[(i - 1) % length, (j - 1) % length]
        l_down = cell[(i + 1) % length, (j - 1) % length]
        r_up = cell[(i - 1) % length, (j + 1) % length]
        r_down = cell[(i + 1) % length, (j + 1) % length]

        neighbours = np.array([up, down, left, right, l_up, l_down, r_up, r_down])
        return neighbours

    def vital_signs(self, cell, i, j):
        neighbours = self.neighbours(cell, i, j)
        alive = 0
        dead = 0

        # find the number of alive/dead neighbors of a given cell
        for n in range(len(neighbours)):
            if neighbours[n] == 1.0:
                alive += 1

            elif neighbours[n] == -1.0:
                dead += 1

        return alive, dead

    def game_of_life(self, cell, i, j):
        """
        Function to apply the Game of Life rules to a target cell
        :param cell: Cell grid
        :param i: x coord of cell
        :param j: y coord of cell
        :return: The value of the cell after the application of the rules
        """

        target_cell = cell[i, j]

        alive, dead = self.vital_signs(cell, i, j)

        if target_cell == 1.0 and alive < 2:
            change = True

        elif target_cell == 1.0 and alive > 3:
            change = True

        elif target_cell == -1.0 and alive == 3:
            change = True

        else:
            change = False

        if change:
            target_cell *= -1

        return target_cell

    def step_forward(self, old_cell):
        """
        Function to apply the Game of Life rules to every cell.

        :param old_cell: The cell grid before the rules are applied.
        :return new_cell: The cell grid after the rules are applied.
        """

        size = self.size
        new_cell = np.ones(np.shape(old_cell))

        for i in range(size):
            for j in range(size):
                target_cell = self.game_of_life(old_cell, i, j)
                new_cell[i, j] = target_cell

        return new_cell

    def get_active_sites(self, cell):
        """
        :param cell: Cell Grid
        :return: The number of alive cells
        """
        size = self.size
        active = 0

        for i in range(size):
            for j in range(size):
                if cell[i, j] == 1:
                    active += 1

        return active

    def find_com(self, cell):
        """
        :param cell: Cell Grid
        :return:
        """
        size = self.size
        n = self.get_active_sites(cell)

        # Change -1 cells to 0
        cell = np.where(cell == 1, cell, cell + 1)
        com_x = 0
        com_y = 0

        for i in range(size):
            for j in range(size):
                com_x += (i * cell[i, j]) / n

        for i in range(size):
            for j in range(size):
                com_y += (j * cell[i, j]) / n

        return com_x, com_y

    def animate_find_stable(self, cell):
        # A test function to visualise finding the stable point
        repeater = 0

        for n in range(self.nstep):
            old_active = self.get_active_sites(cell)

            cell = self.step_forward(cell)
            new_active = self.get_active_sites(cell)

            if old_active == new_active:
                repeater += 1
            else:
                repeater = 0

            if repeater == 10:
                print('IT IS STABLE!!!')
                print('Stabilised after ' + str(n) + ' iterations.')

                f = open('equilibration_time.txt', 'a')
                f.write(str(n) + ',\n')
                f.close()
                break
            #       show animation
            plt.cla()
            im = plt.imshow(cell, cmap='Greens', animated=True)

            plt.draw()

            plt.pause(0.0001)

    def short_find_stable(self, cell):
        """
        A function to find the point of equilibration
        :param cell: cell grid
        """
        repeater = 0

        for n in range(self.nstep):
            old_active = self.get_active_sites(cell)

            cell = self.step_forward(cell)
            new_active = self.get_active_sites(cell)

            if old_active == new_active:
                repeater += 1
            else:
                repeater = 0

            if repeater == 10:
                print('Stabilised after ' + str(n) + ' iterations.')
                f = open('equilibration_time.txt', 'a')
                f.write(str(n) + ',\n')
                f.close()
                break

    def get_equilibrium_times(self, iterations):
        for i in range(iterations):
            print(i)
            cell = self.init_cells_random()
            self.short_find_stable(cell)

    def com_data(self, cell):
        """
        Function to find the coordinates of the centre of mass and write to a file
        :param cell: cell grid
        """
        f = open('com.txt', 'w')
        bc = self.size / 2

        for n in range(self.nstep):
            cell = self.step_forward(cell)

            glider_pos = np.where(cell == 1)

            # checking bc for x
            glider_x = glider_pos[0]
            x_min = np.min(glider_x)
            x_max = np.max(glider_x)
            x_test = x_max - x_min

            # checking bc for y:
            glider_y = glider_pos[1]
            y_min = np.min(glider_y)
            y_max = np.max(glider_y)
            y_test = y_max - y_min

            if x_test < bc and y_test < bc:
                com_x, com_y = self.find_com(cell)
                f.write('%i %f %f\n' % (n, com_x, com_y))
            else:
                continue

        f.close()

    def animate(self, cell):
        """
        Function show an animation of the system
        :param cell: The cell grid
        """

        for n in range(self.nstep):
            cell = self.step_forward(cell)

            #       show animation
            plt.cla()
            im = plt.imshow(cell, cmap='Greens', animated=True)

            plt.draw()

            plt.pause(0.0001)

'''
Uncomment to calculate quantities and write to file.
'''
# gol.animate_find_stable(random_cells)
# gol.com_data(glider_cells)
# gol.get_equilibrium_times(10)
