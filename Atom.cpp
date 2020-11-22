#include "Atom.h"



using namespace std;

array<int,2> Atom::get_r(){
	return r;
}

//Functions definition of Atom class:

void Atom::move(array<double,6> probaMove,array<array<int,2>,6> indexNeighbours, double randomNumber){
	
	double N1=probaMove[0]; //Put together the probabilities to draw only one random number per move
	double N2=N1+probaMove[1];
	double N3=N2+probaMove[2];
	double N4=N3+probaMove[3];
	double N5=N4+probaMove[4];
	double N6=N5+probaMove[5];
	if (N6>1.0){ //Check if probability> => normalization + cout message
		cout<<N6<<endl;
		N1/=N6;
		N2/=N6;
		N3/=N6;
		N4/=N6;
		N5/=N6;
		N6=1.0;
		cout<<"Total mouvement probability over 1, renormaization happened"<<endl;
	}
	if(randomNumber<N1){
		isMoving=true;
		r=indexNeighbours[0];
		}
	else if(randomNumber<N2){
		isMoving=true;
		r=indexNeighbours[1];
		}
	else if(randomNumber<N3){
		isMoving=true;
		r=indexNeighbours[2];
		}
	else if(randomNumber<N4){
		isMoving=true;
		r=indexNeighbours[3];
		}
	else if(randomNumber<N5){
		isMoving=true;
		r=indexNeighbours[4];
		}
	else if(randomNumber<N6){
		isMoving=true;
		r=indexNeighbours[5];
		}

	else{isMoving=false;}

}
bool Atom::isMovingStatus(){
	if(isMoving==true){return true;}
	else if(isMoving==false){return false;}
	else{ cout<<"isMoving not declared, can't assign true or false value"<<endl;
	return false;
	}

}
//bool Atom::isTouchingStatus(){return isTouching;}
