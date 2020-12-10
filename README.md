# dendriticGrowth
Simulation project: Symmetries in dendritic growth on a Platinum (111) substrate.
Marc Jacquart

This program has been created on ubuntu (linux), libraries and program needed:
c++17(for directory creation and random library), g++, python3, matplotlib, numpy

READ ME:

To exectute: Use makefile to create to compile the program (1) and launch it (2) with comands:
(1) make
(2) ./main

This program uses a kinetic Monte Carlo algorithm to compute the evolution
of a atomic system in order to reproduce epiaxial growth of atomic clusters.

Change global physical parameters in main() function in main.cpp
Change energy barriers in System.h class
Change allowed movements in computeProba() function in System.cpp

By setting the boleean setting in addAtom() to false, one can decide that a newly generated atom 
will not stick to other (non stick) atoms. If this atom sticks comes beside a "sticky" atom,
it becomes sticky itself. This allows the user to chose the number and position of cluster by placing
a few sticky atoms at the beginning of the simulation.
This option is unphysical but handy to study the shape of the cluster with far less iterations by 
imposing only one final custer.

Results are found in test files in Results/Current folder. 
Change the directory name before running a new simulation to avoid overwriting.
This folder contains 3 text files: 
	Parameter.txt: 	Contains a summary of the used parameters during the simulation
	Simulation.txt:	Contains the results of the simulation: Each line correspond
					to the position of the atom at a given timestep (to increase performences, 
					the the print frequency can be set in main.cpp)
	Time.txt:		Contains the time between each timestep for performance analysis.


Analysis is done by the Analysis.py script. 
Uncomment the function calls at the end of the program to chose which analysis to do:
	plotGif() 		return a gif animation with all the timesteps of Simulation.txt
	plotLast() 		return a plot of the last frame
	plotDistance() 	return an histogram of all the distances between the atoms 
					(try to find cluster distances and sizes)
 


Indices of surrounding fcc sites around an atom A (Used in all the tab in System.cpp):

   6    1
	  \  /
   	 \/
5----A----2
  	 /\
  	/  \
   4    3
