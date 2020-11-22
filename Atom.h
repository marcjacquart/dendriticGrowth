#pragma once
#include <iostream>
#include <array>
using namespace std;

class Atom{
	public:
	
	Atom(array<int,2> position, bool isStk=false ) //Constructor
		:r(position),  isSticking(isStk) ,isMoving(true)//isTouching(false),
		{/*isSticking=isStk;*/
		//cout<<"Atom cerated at X: "<<position[0]<<" et Y: "<<position[1]<<" and isStucking: "<<isStk<<endl;
	}
	
	void update(array<double,2>newPos);

	array<int,2> get_r();
	bool isMovingStatus();
	//bool isTouchingStatus();

	void move(array<double,6> probaMove,array<array<int,2>,6> indexNeighbours, double randomNumber); //return new position
	bool getSticking(){//return true;
		if(isSticking==true){/*cout<<"TRUE"<<endl;*/return true;}
		else if (isSticking==false){/*cout<<"FALSE"<<endl;*/return false;}
		else{cout<<"ERROR GET STICKING at X:"<<get_r()[0]<<" Y: "<<get_r()[1]<<endl;
			return false;}
		}
	void setSticking(bool valueStick){isSticking=valueStick;/*cout<<"Stick updated"<<endl;*/}

	private:

	array<int,2> r;			//int car coordonnees

	bool isSticking;
	bool isMoving;
	//bool isTouching;
	//array<double,6>probaMvmt;

	};
