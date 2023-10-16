####################################################################################################

To show the animation, run the Main.py file then follow the input prompts.

- For Game of Life, input "G [grid size] [n steps]". Tou will then be prompted to input 1, 2 or, 3 to signify whether
  the initial state should be random, a blinker, or a glider.

- For SIRS, input "S [grid size] [n steps] [p1] [p2] [p3] [immune fraction]"
- If you wish to display the animation without an immune population, set the immune fraction to 0.

- Probabilities for different SIRS states:
    Absorbing state: 0.1, 0.5, 0.1
    Dynamic Equilibrium: 0.5, 0.5, 0.5
    Waves: 0.9, 0.6, 0.1

####################################################################################################

Calculations on the Systems are done by uncommenting the relevant functions at the bottom of their files and running
the individual files, GameOfLife.py or SIRS.py

Graphs are produced by the functions in the Graphs.py file.

NOTE: To save time, there are pre-produced graphs that can be found in the Plots folder. These are produced by the code
as it in the current version.