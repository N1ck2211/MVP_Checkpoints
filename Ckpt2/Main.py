from SIRS import SIRS
from GameOfLife import GameOfLife


def main():
    print()
    print('If you wish to run the Game of Life please input G, followed by integer values of grid size and no. steps.')
    print('If you wish to run the SIRS model please input S, followed by values of grid size, no. steps, p1, '
          'p2, p3, and the fraction of immune cells')

    print()

    task = input('Input task (G/S) and Values: ')

    task = task.split(' ')

    model = task[0]
    print(task[1])

    # param_three = task[3]
    if model == 'G' or model == 'g':
        print('You have selected the Game of Life')

        size = int(task[1])
        nstep = int(task[2])

        setup = input('Random, blinker, or glider setup? (1/2/3) ')

        if setup == '1':
            gol = GameOfLife(size, nstep)
            random_cells = gol.init_cells_random()
            gol.animate(random_cells)

        elif setup == '2':
            gol = GameOfLife(size, nstep)
            blinker = gol.init_cells_blinker()
            gol.animate(blinker)

        elif setup == '3':
            gol = GameOfLife(size, nstep)
            glider = gol.init_cells_glider()
            gol.animate(glider)

        else:
            print('Invalid input, sorry')

    elif model == 'S' or model == 's':
        print('You have selected SIRS')
        size = int(task[1])
        nstep = int(task[2])
        p1 = float(task[3])
        p2 = float(task[4])
        p3 = float(task[5])
        im_frac = float(task[6])

        sirs = SIRS(size, nstep, p1, p2, p3, im_frac)
        cells = sirs.init_cells_immune_frac()

        sirs.animate_immune(cells)

    else:
        print('Invalid input, sorry')


main()

'''
Absorbing state: 0.1, 0.5, 0.1
Dynamic Equilibrium: 1, 1, 1
Waves: 0.9, 0.6, 0.1
'''
