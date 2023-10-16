Note: if the code is run to graph,

#####################################################################################################################

- To run the file to create an animation or produce data, run main.py with the run() function uncommented.
- You will be prompted first whether you want to produce an animation or a data plot. Press 'A' to produce an animation
  or press 'G' to produce a plot.

- If you choose to create a plot:
- You will then be prompted to input the dynamics used. Press 'G' to use glauber dynamics, press 'K' to use Kawasaki.
- You then input the length of the matrix, N, which will be used to produce the NxN matrix.
- Then input the step size of the temperature, this will be the temperature increments between 1 and 3.
- Finally, choose the desired quantity to be graphed, input 'M' for Magnetisation, 'S' for susceptibility, 'E' for
  energy, and 'C' for heat capacity, followed by the number of steps the program should run for.

- If you choose to animate:
- First you will be prompted to input the size of the matrix, the temperature of teh system, and the number of steps
  you would like to run for.
- You will then be prompted to input the dynamics used. Press 'G' to use glauber dynamics, press 'K' to use Kawasaki.

######################################################################################################################

- If you wish to plot graphs from the data present in the output files, comment out the run() function at the end
 of the main file, and uncomment the block of functions designed to plot from the various .txt files.