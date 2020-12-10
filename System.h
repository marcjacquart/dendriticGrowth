#pragma once
#include <iostream>
#include "Atom.h"
#include <array>
#include <vector>
#include <fstream>
#include <random>
using namespace std;
class System{
	
	public:
	System(double Temperature, int latticeSizeX, int latticeSizeY, int maxAtomsInit, mt19937 randomGenerator)
	:nAtomsInSyst(0),												// Fill all the needed parameters in the initialization list
	T(Temperature),
	sizeX(latticeSizeX),
	sizeY(latticeSizeY),
	atomList(maxAtomsInit,Atom({-2,-2})),
	maxAtoms(maxAtomsInit),
	generator(randomGenerator), 
	fcc(latticeSizeX,vector<Atom*>(latticeSizeY,NULL)){				// The whole tab is filled with 0
			
		double eToProba=energyToProba(E_M,1.0); 					// Need to set dt=1.0 for energyToProba in order to get the real dt
		cout<<"ENERGYTOPROBA E_M (dt=1): "<<eToProba<<endl;
		dt=fMoveStep/(6*constPropExp*exp(-E_M/(k_b*T)));			// Compute optimized dt to allow a free movement frequency of fMoveStep in average
		cout<<"dt= "<<dt<<endl;	
		cout<<"ENERGYTOPROBA E_M(bon dt): "<<energyToProba(E_M,dt)<<endl;
		uniform_real_distribution<> distribMove(0.0,1.0);
		nStepsTot=100/dt;											// Change 100 if time simulated is not 100s. Just for information, won't take the time from main.cpp for simplicity
		cout<<"System created, will be updated every dt="<<dt<<" s for a maximum of 1/10 atoms moving at each step."<<endl;
		cout<<"This means the simulation will contain"<<nStepsTot<<" steps."<<endl;

		energyToProbaTab[0]=energyToProba(E_M,dt);					// Fill the proba tab with all the probability calculated from the energy barrier for faster use in the system evolution
		energyToProbaTab[1]=energyToProba(E_Ac,dt);
		energyToProbaTab[2]=energyToProba(E_Bc,dt);
		energyToProbaTab[3]=energyToProba(E_Ae,dt);
		energyToProbaTab[4]=energyToProba(E_Be,dt);
		energyToProbaTab[5]=energyToProba(E_Aec,dt);
		energyToProbaTab[6]=energyToProba(E_Bec,dt);
		for(int j(0);j<7;j++){cout<<"Proba "<<j<<" = "<<energyToProbaTab[j]<<endl;}
		


	}
	

	// Method and attributes definition
	array<double,6> computeProba(Atom const &A);						// Compute the probability movement of an atom for the 6 surounding positions
	bool addAtom(array<int,2> const (&pos), bool isSticking=false);		// Add an atom to the system, not sticking by default
	void printSystem(string const &outFile, bool printAdatom=true);		// Print the simulation results: position of the atoms in the system
	array<array<int,2>,6> indexNeighbours(array<int,2> const(&pos));	// index {X,Y} of each neighbour of one given Atom
	void moveSystem();													// Move every atom in the atomList
	double get_dt(){return dt;}
	void printParam(string const &paramFile);							// Print physical parameters in the text file at the begining
	int get_nStepsTot(){return nStepsTot;}
	private:
		
	int nAtomsInSyst;
	void printFccStatus();							// For debug: print all the position of all atoms in the fcc pointer tab (to verify if they are not corrupted)
	//Attributes:
	array<double,7>energyToProbaTab;
	double const fMoveStep=1.0; 					// Fraction of total atoms that move at each timestep
	double const T;
	int const sizeX;
	int const sizeY; 								// Lattice size
	int nStepsTot;
	vector<Atom> atomList; 							// List of all atoms in the system
	int const maxAtoms;

	//Random:
	mt19937 generator;
	uniform_real_distribution<> distribMove;

	vector<vector<Atom*>>fcc; 						// Efficient way to keep trace of the atoms(propability calcul with neighbours). Only on fcc sites are considered in this model
	double const constPropExp=pow(10,12); 			// Proportionnality constant, 1.0 for now,  will be corelated to dt
	
	double const k_b=8.617e-5;						// Boltzmann constant [eV/K]

	//Energy barriers: Calculated by effective medium theory
	double const E_M=0.08;							//[eV] Energy barrier for normal fcc displacement
	double const E_Ac=0.08;							//[eV] Energy barrier for corner A step
	double const E_Bc=0.248;						//[eV] Energy barrier for corner B step
	double const E_Ae=0.187;						//[eV] Energy barrier for edge A step
	double const E_Be=0.389;						//[eV] Energy barrier for edge B step
	double const E_Aec=0.244;						//[eV] Energy barrier for edge to corner A step
	double const E_Bec=0.415;						//[eV] Energy barrier for edge to corner B step
 	

	//Methods:
	int numberOfNeighbours (array<int,2> const (&pos)); 										// Check if there is something in each of the 6 surounding spots to return the  number of neighbouring atoms of one site
	void updateAtomPointer(Atom &A, array<int,2>const (&oldPos), array<int,2> const (&newPos)); // Update the fcc pointer 2d array
	double energyToProba(double energyBarrier,double dt); 										// Used at the begining to compute the probability of movement from the energy barriers
	double dt;
	array<bool,6> posNeighbours(array<int,2> const (&pos)); 									// Index 1 to 6 with respect to convention (see README.txt), true if an atom is here
	void updateStick(Atom &A); 																	// Set isSticking to true if the atom is beside of that is sticking
	
	
};
