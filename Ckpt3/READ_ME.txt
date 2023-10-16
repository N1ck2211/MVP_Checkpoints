######################################################################################################

- To run the program run the main.py file. You will then be prompted to input the algorithm you want to run.
- To run the Cahn-Hilliard, input C
- To run the poisson equation for an electric field press E
- To run the poisson equation for a magnetic field due to a wire press M
- To run the Gauss-Seidel equation for an electric field press G
- To see the convergence of the SOR algorith press S

######################################################################################################

- When running hte Cahn-Hilliard, you will be prompted to input whether you wish to view an animation, or gather the
 data for a plot of the free energy. For both you will be prompted to input the array size, the number of iterations,
 and the value of phi_o. Phi_o = 0 will give spinodal decomposition, while 0.5 gives droplets. A recommended step count
 for the free energy graph is 100000.

######################################################################################################
-For the field calculations, you will be asked to input the size of the array, the max number of iterations, and the
tolerance. A recommended tolerance is 0.01.

- When tolerance is reached, a contour plot of the potential is displayed, as well as a vector plot of the given field.

- Data files for the magnitude of the potential and field against distance to the midpoint are saved automatically.

- When the SOR is run, a datafile of the omega vs time to converge is produced.

######################################################################################################

- To view the graphs of the free energy, potentials, fields, and sor convergence, run the GraphFromFile.py file.
- For convenience, the outputs of this file are displayed in the Plots directory.