#pragma once
#include <iostream>
#include <array>
using namespace std;

class Atom{
	public:
	
	Atom(array<int,2> const (&position), bool isStk=false ) 						// Constructor
		:r(position),  isSticking(isStk) ,isMoving(true) 	// isTouching(false), 	// For cluster individual size analysis
		{}
	


	array<int,2> get_r()const;
	bool isMovingStatus()const; 			// Getters don't modify the object -> const functions
	//bool isTouchingStatus();

	void move(array<double,6> const (&probaMove),array<array<int,2>,6> const (&indexNeighbours), double randomNumber); //return new position

	bool getSticking() const{				// const because don't change the object
		if(isSticking==true){return true;}
		else if (isSticking==false){return false;}
		else{cout<<"ERROR GET STICKING at X:"<<get_r()[0]<<" Y: "<<get_r()[1]<<endl;
			return false;}
		}
		
	void setSticking(bool valueStick){isSticking=valueStick;}

	private:

	array<int,2> r;							// Coordinates are integers, need an array of 2 int

	bool isSticking;
	bool isMoving;
	//bool isTouching;						// Could be useful to separate the cluster for individual size analysis, not done at the moment
	//array<double,6>probaMvmt;

	};
