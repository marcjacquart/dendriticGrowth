/* Simulation Project:
 * 
 * Symetries in dendritic 
 * growth on a Platinum 
 * (111) substrate.
 * 
 * Marc Jacquart
 */

#include <iostream>
#include "System.h"
#include <random>
#include <stdio.h> //To remove the txt file at the beginning

using namespace std;



int main(){
	//Parameters Init: 
	double const T(36.0); 				//[K]  Temperature
	//double const dt(0.1); 			//[s]  Timestep
	//double const tFinalDrop(1); 		//[s]  After this no more atoms are generated
	double const tFinal(100);			//[s]  End of the simulation
	//double const nAtoms(6000);			//[-]  Total number of diffusing atoms
	int const latticeSizeX(400);		//[-]  Size of the lattice grid in X
	int const latticeSizeY(400);		//[-]  Size of the lattice grid in y
	double totalCoverage(0.04);			//Fraction of the lattice we want covered with deposited atoms
	int printFrequency=100000;			//Will print all timesteps in the txt file if true, only the last if false
	//Energy barriers: directly set in System.h for easier declaration of a system

	string const outputFile("Simulation.txt");
	remove(outputFile.c_str());
	cout<<"Symetries in dendritic  growth on a Platinum (111) substrate. Marc Jacquart"<<endl;
	cout<<endl;
	cout<<"Time simulated: "<<tFinal<<" [s]"<<endl;
	cout<<"Parameters:"<<endl;
	cout<<"Temperature: "<<T<<" [K]"<<endl;
	cout<<"Fraction of plane filled with atoms: "<<totalCoverage<<endl;
	cout<<"lattice size X: "<<latticeSizeX<<" and Y: "<<latticeSizeX<<endl;
	
	//Random engine: from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution

	random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<> distribAdd(0.0,1.0);
	uniform_int_distribution<> distribCoordX(0,latticeSizeX-1);
	uniform_int_distribution<> distribCoordY(0,latticeSizeY-1); //6 significant digits
	

	//Part1: System init
	vector<Atom> atomTab;		// Atom.h is called by system.h Create empty tab for the begining
	System S(T,latticeSizeX,latticeSizeY,atomTab, generator);			// System is created (all diffusing atoms)
	//Add the first sticking thetraedre
	int coordX =latticeSizeX/3.5;
	int coordY =latticeSizeY/2;
	//cout<<coordX<<" "<<coordY<<endl;
	S.addAtom({coordX,coordY},true);
	S.addAtom({coordX+1,coordY},true);
	S.addAtom({coordX,coordY+1},true);
	S.addAtom({coordX+1,coordY+1},true);
	//S.addAtom({0,0},true);

	//Calculate number of steps to simulate:
	long int nSteps=tFinal/S.get_dt();
	cout<<"Number of steps: "<<nSteps<<"with dt= "<<S.get_dt()<<"and tFinal="<<endl;
	//probability of atom deposition:
	long double propaDeposit = (totalCoverage * latticeSizeX * latticeSizeY)/nSteps;
	cout<<"Probability of atom deposition at each timestep: "<<propaDeposit<<endl;
	//Part 2: Deposit the atoms

	for(int step(0);step<nSteps;step++){ //nAtomSyst: Number of atomas in the system: 0 at the beginning
		double N= distribAdd(generator); //number of atoms to put this timestep: 0 or 1
		if (N<propaDeposit){
			int X= distribCoordX(generator);

			int Y= distribCoordY(generator);
			//if(S.addAtom({X,Y}, false)==true){nAtomSyst++;}; //many clusters if true here //To have a fix number of final atoms
			S.addAtom({X,Y}, false); //Adding an atom that wont stick at first
			}
		S.moveSystem();
		if (step%printFrequency==0){S.printSystem(outputFile);cout<<step<<endl;} //Same configuration writting can be optimized
		}
	/*cout<<"All atoms have been deposited"<<endl;
	for (int j(0);j<100000;j++){
		if(j%10000==0){cout<<j<<endl;}
		S.moveSystem();
		if (printAll==true){S.printSystem(outputFile);}
	}*/ //Plus besoin de bouger sans deposer
	S.printSystem(outputFile); //will anyway print the last configuration
	cout<<"End of the simulation"<<endl;
	return 0;
}
