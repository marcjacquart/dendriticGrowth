#include <iostream>
#include "System.h"
#include <random>						// Generate random number
#include <chrono> 						// Elapsed time
#include <ctime>						// To generate directory name for the new simulation
#include <string>
#include <filesystem> 					// Create directory
#include <stdlib.h>						// For exit failure
using namespace std;


int main(){
	//Initial Parameters: 
	double const T(30.0); 				// [K]  Temperature
	double const tFinal(100);			// [s]  End of the simulation
	int const latticeSizeX(300);		// [-]  Size of the lattice grid in X
	int const latticeSizeY(300);		// [-]  Size of the lattice grid in y
	double totalCoverage(0.06);			// Fraction of the lattice we want covered with deposited atoms
	bool printAdatom(true);				// Do we want to print the non sticking
	int printFrequency=1000000;			// Will print all timesteps in the txt file if true, only the last if false
	//Energy barriers: directly set in System.h for easier declaration of a system

	//double const nAtoms(6000);		// [-]  Total number of diffusing atoms //These parameters are calculated directly
	//double const dt(0.1); 			// [s]  Timestep
	//double const tFinalDrop(1); 		// [s]  After this no more atoms are generated

	//Starting chrono:
	cout<<"Starting chrono..."<<endl;
	auto start= chrono::system_clock::now();


	//Output test files
	time_t dirTime = chrono::system_clock::to_time_t(start);

	//change time format:
	struct tm * timeinfo;
  	char buffer [80];
  	timeinfo = localtime (&dirTime);

  	strftime (buffer,80,"%d_%m_%G_%R_",timeinfo);							// Setting time format
  	string formattedTime = buffer;
  	double totalCoverageName = 100*totalCoverage;							// Coverage in percent
  	string const directoryName("./Results/"+formattedTime+to_string(latticeSizeX)+"x"+to_string(latticeSizeY)+"_T"+to_string((int)T)+"_"+to_string((int)totalCoverageName)+"%");

  	if(filesystem::create_directories(directoryName)==false){cout<<"Error while creating new directory for the simulation"<<endl;exit(EXIT_FAILURE);} // Can happend if two simulations happens to have the same name. Better to exit at the begining
	cout<<"Internal Clock time: "<<formattedTime<<endl;

	string const outputFile(directoryName+"/Simulation.txt");				// Creating the 3 output files
	string const outputTime(directoryName+"/Time.txt");
	string const outputParam(directoryName+"/Parameters.txt");



	cout<<"Symetries in dendritic  growth on a Platinum (111) substrate. Marc Jacquart"<<endl; // Cout the parameters at the begining of the simulation
	cout<<endl;
	cout<<"Time simulated: "<<tFinal<<" [s]"<<endl;
	cout<<"Parameters:"<<endl;
	cout<<"Temperature: "<<T<<" [K]"<<endl;
	cout<<"Fraction of plane filled with atoms: "<<totalCoverage<<endl;
	cout<<"lattice size X: "<<latticeSizeX<<" and Y: "<<latticeSizeX<<endl;

																								// Copy these infos in the parameter text file for further use:
	ofstream writeParam;																		// Define output stream
	writeParam.open(outputParam,ios_base::app); 												// Documment with parameters
	writeParam<<"Parameters:"<<endl<<endl;
	writeParam<<"Time simulated: "<<tFinal<<" [s]"<<endl;
	writeParam<<"Temperature: "<<T<<" [K]"<<endl;
	writeParam<<"Fraction of plane filled with atoms: "<<totalCoverage<<endl;
	writeParam<<"lattice size X: "<<latticeSizeX<<" and Y: "<<latticeSizeX<<endl;
	writeParam<<"Print Frequency: "<<printFrequency<<endl;
	writeParam.close();
	
	//Random engine: from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution
	random_device rd;  																			// Will be used to obtain a seed for the random number engine
    mt19937 generator(rd()); 																	// Standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<> distribAdd(0.0,1.0);											// Create 6 significant digits pseudo random number between [0.0,1.0]
	uniform_int_distribution<> distribCoordX(0,latticeSizeX-1);
	uniform_int_distribution<> distribCoordY(0,latticeSizeY-1); 								
	

	//Part1: System init
	
	int maxAtoms (latticeSizeX*latticeSizeY*totalCoverage*1.01);								// Max allowed atoms in he system (add 1% for deposit fluctuations) to declare atomList as an array and stop this pointer nonesense when pushing back :)
	System S(T,latticeSizeX,latticeSizeY,maxAtoms/*,atomTab*/, generator);						// System is created (all diffusing atoms)
	S.printParam(outputParam); 																	//Print System parameters in txt file
	//Add the first sticking thetraedre for one cluster simulation, then S.addAtom({X,Y}, false); in the loop
	/*
	int coordX =latticeSizeX/3.5;
	int coordY =latticeSizeY/2;
	S.addAtom({coordX,coordY},true); //bug 0,0
	S.addAtom({coordX+1,coordY},true);//bug
	S.addAtom({coordX,coordY+1},true);
	S.addAtom({coordX+1,coordY+1},true);
	*/

	//Calculate number of steps to simulate:
	long int nSteps=tFinal/S.get_dt();
	cout<<"Number of steps: "<<nSteps<<"with dt= "<<S.get_dt()<<"and tFinal="<<endl;
	//probability of atom deposition:
	long double propaDeposit = (totalCoverage * latticeSizeX * latticeSizeY)/nSteps;
	cout<<"Probability of atom deposition at each timestep: "<<propaDeposit<<endl;
	
	//Part 2: Deposit the atoms
	for(int step(0);step<nSteps;step++){ 														// For all the time steps of the simulation:
		double N= distribAdd(generator); 														//	Random number generation to tell if an atom is aded to the system
		if (N<propaDeposit){
			int X= distribCoordX(generator);
			int Y= distribCoordY(generator);
			
			S.addAtom({X,Y}, true); 															// Adding an atom that (true=will, false =wont) stick at first
			}
		S.moveSystem();
		if (step%printFrequency==0){
			S.printSystem(outputFile);cout<<step<<endl;
			auto end= chrono::system_clock::now();
			ofstream writeTime;																	// Define output Time file stream
			writeTime.open(outputTime,ios_base::app); 											// Text documment with results
			cout<<"Elapsed time: "<<chrono::duration_cast<chrono::seconds>(end-start).count()<<" [s]."<<endl;
			double percentageDone=100*( (double)step/S.get_nStepsTot())*( (double)step/S.get_nStepsTot());
			cout<<"Simulation progress: "<<(int)percentageDone<<" %"<<endl;						// Given the linear dependency of time to number of atoms in the system. Integer for display
			writeTime<<chrono::duration_cast<chrono::seconds>(end-start).count()<<endl;
			writeTime.close();
			} 
		}
	
	S.printSystem(outputFile, printAdatom); 													// Will anyway print the last configuration
	cout<<"End of the simulation"<<endl;
	return 0;																					// End of the simulation
}
