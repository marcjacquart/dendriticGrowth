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
	System(double Temperature, int latticeSizeX, int latticeSizeY, vector<Atom>listOfAtoms, mt19937 randomGenerator)
	:T(Temperature),
	sizeX(latticeSizeX),
	sizeY(latticeSizeY),
	atomList(listOfAtoms),
	generator(randomGenerator), 
	fcc(latticeSizeX,vector<Atom*>(latticeSizeY,NULL))
	
	{
	cout<<"dt= "<<dt<<endl;							//Will change later, needed 1.0 for energyToProba in order to get dt
	double eToProba=energyToProba(E_M,1.0);
	cout<<"ENERGYTOPROBA E_M (dt=1): "<<eToProba<<endl;
	dt=fMoveStep/(6*constPropExp*exp(-E_M/(k_b*T)));
	cout<<"dt= "<<dt<<endl;	
	cout<<"ENERGYTOPROBA E_M(bon dt): "<<energyToProba(E_M,dt)<<endl;
	uniform_real_distribution<> distribMove(0.0,1.0);
	
	cout<<"System created, will be updated every dt="<<dt<<" s for a maximum of 1/10 atoms moving at each step."<<endl;
	cout<<"This means the simulation will contain"<<100/dt<<" steps."<<endl;
	}
	


	array<double,6> computeProba(Atom &A);
	bool addAtom(array<int,2>pos, bool isSticking=false);
	void printSystem(string outFile);
	array<array<int,2>,6> indexNeighbours(array<int,2>pos);
	void moveSystem();
	double get_dt(){return dt;}

	private:


	//Attributes:
	double const fMoveStep=0.99; //fraction of total atoms that move at each timestep
	double const T;
	int sizeX;
	int sizeY; 							//lattice size
	vector<Atom> atomList; 				//list of all atoms in the system

	//Random:
	mt19937 generator;
	uniform_real_distribution<> distribMove;

	vector<vector<Atom*>>fcc; 			//efficient way to keep trace of the atoms(propability calcul with neighbours). For now only on fcc sites
	//vector<vector<Atom*>>fcc;
	double const constPropExp=pow(10,12); 		//proportionnality constant, 1.0 for now,  will be corelated to dt
	
	double const k_b=8.617e-5;			//Boltzmann constant [eV/K]

	//Energy barriers: Calculated by effective medium theory
	double const E_M=0.08;				//[eV] Energy barrier for normal fcc displacement
	double const E_Ac=0.08;				//[eV] Energy barrier for corner A step
	double const E_Bc=0.248;			//[eV] Energy barrier for corner B step
	double const E_Ae=0.187;			//[eV] Energy barrier for edge A step
	double const E_Be=0.389;			//[eV] Energy barrier for edge B step
	double const E_Aec=0.244;			//[eV] Energy barrier for edge to corner A step
	double const E_Bec=0.415;			//[eV] Energy barrier for edge to corner B step
 	

	//Methods:
	int numberOfNeighbours (array<int,2>pos);
	void updateAtomPointer(Atom& A, array<int,2>oldPos,array<int,2>newPos); //update the pointer grid
	double energyToProba(double energyBarrier,double dt);
	double dt;
	array<bool,6> posNeighbours(array<int,2>pos);
	void updateStick(Atom &A);
	
	
};
