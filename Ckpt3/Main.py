from ZweiFisch import Poisson
from MagneticFisch import MagneticPoisson
from GaussSeidel import GaussSeidel
from SOR import SOR, convergence_test
from CahnHilliard import CahnHillard


def run():
    algo = input('Which algorith would you like to run? (C/E/M/G/S) ')

    if algo == 'C' or algo == 'c':
        decision = input('Would you like to see an animation, or a plot of the free energy? (A/P) ')

        if decision == 'A' or decision == 'a':
            inputs = input('Input [array size] [step no] [phi_0]: ')
            inputs = inputs.split(' ')
            size = int(inputs[0])
            nstep = int(inputs[1])
            phi_o = float(inputs[2])
            ch = CahnHillard(size, nstep, phi_o)

            ch.initialise_system()
            ch.animate()

        if decision == 'P' or decision == 'p':
            inputs = input('Input [array size] [step no] [phi_0]: ')
            inputs = inputs.split(' ')
            size = int(inputs[0])
            nstep = int(inputs[1])
            phi_o = float(inputs[2])

            ch = CahnHillard(size, nstep, phi_o)

            ch.initialise_system()
            ch.get_f()

        else:
            print('Invalid Input')

    if algo == 'E' or algo == 'e':
        inputs = input('Input [array size] [step no] [tolerance]: ')
        inputs = inputs.split(' ')
        size = int(inputs[0])
        nstep = int(inputs[1])
        tolerance = float(inputs[2])

        psn = Poisson(size, nstep, tolerance)
        psn.initialise()
        # psn.animate()
        psn.endpoint()

    if algo == 'M' or algo == 'm':
        inputs = input('Input [array size] [step no] [tolerance]: ')
        inputs = inputs.split(' ')
        size = int(inputs[0])
        nstep = int(inputs[1])
        tolerance = float(inputs[2])

        mp = MagneticPoisson(size, nstep, tolerance)
        mp.initialise()
        mp.endpoint()

    if algo == 'G' or algo == 'g':
        inputs = input('Input [array size] [step no] [tolerance]: ')
        inputs = inputs.split(' ')
        size = int(inputs[0])
        nstep = int(inputs[1])
        tolerance = float(inputs[2])

        gs = GaussSeidel(size, nstep, tolerance)
        gs.initialise()
        # gs.animate()
        gs.endpoint()

    if algo == 'S' or algo == 's':
        # sor = SOR(50, 100000, 1.5, 0.01)
        # sor.initialise()
        # sor.endpoint()
        convergence_test()

    else:
        print('Invalid Input')


run()
