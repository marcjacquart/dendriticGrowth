#include "System.h"
#include<cmath>
using namespace std;
//Methods for the System Class:

void System::updateAtomPointer(Atom& A, array<int,2>oldPos,array<int,2>newPos){
	
	if(oldPos[0]!=-1){fcc[oldPos[0]][oldPos[1]]=NULL;/*cout<<"old pointeer destroyed"<<endl;*/}
	
	fcc[newPos[0]][newPos[1]]=&A; //assign in tab the adress of A, wich is a reference (just another name for the variable) over the atom we want
	//cout<<"updateAtomPointer for an atom at X: "<<fcc[newPos[0]][newPos[1]]->get_r()[0]<<" and Y: "<<fcc[newPos[0]][newPos[1]]->get_r()[1]<<endl;
}


void System::printSystem(string outFile){
	ofstream writeTxt;				//Define output file stream
	writeTxt.open(outFile,ios_base::app); //Text documment with results
	if (not writeTxt.fail()){
		for(size_t i(0);i<atomList.size();i++){
			double coordX(atomList[i].get_r()[0]);
			double coordY(atomList[i].get_r()[1]);
			writeTxt<< fmod(coordX+0.5*coordY,sizeX)<<","<< fmod(0.866*coordY,sizeY)<<";"; //gives: x_0,y_0;x_1,y_1;x_2,y_2;... fmod:modulo for double

		}
		writeTxt<<endl; //one line for each timestep
	}
	writeTxt.close();
}




array<array<int,2>,6> System::indexNeighbours(array<int,2>pos){
	
	array<int,2>index1={pos[0],(pos[1]+1)%sizeY};
	array<int,2>index2={(pos[0]+1)% sizeX,pos[1]};
	array<int,2>index3={(pos[0]+1)%sizeX,(pos[1]-1+sizeY)%sizeY}; //+sizeY to avoid -1%sizeY=-1 :/ reminder not modulo
	array<int,2>index4={pos[0],(pos[1]-1+sizeY)%sizeY};
	array<int,2>index5={(pos[0]-1+sizeX)%sizeX,pos[1]};
	array<int,2>index6={(pos[0]-1+sizeX)%sizeX,(pos[1]+1)%sizeY};
	return{index1,index2,index3,index4,index5,index6};

}


int System::numberOfNeighbours (array<int,2>pos){
	array<array<int,2>,6> index =indexNeighbours(pos);
	int result(0);
	for(int i(0);i<6;i++){
		//cout<<"INDEX: "<<index[i][0]<<" , "<<index[i][1]<<endl;
		if(fcc[index[i][0]][index[i][1]]!=NULL){result+=1;}
	}
	return result;

}
array<bool,6> System::posNeighbours(array<int,2>pos){ //index 1 to 6 with respect to convention, true if an atom is here
	array<bool,6> results ={false};
	array<array<int,2>,6> index =indexNeighbours(pos);
	for(int i(0);i<6;i++){
		if(fcc[index[i][0]][index[i][1]]!=NULL){results[i]=true;}
	}
	return results;
}

bool System::addAtom(array<int,2>pos, bool isSticking){ //not sticking by default
	if (fcc[pos[0]][pos[1]]!=NULL){return false;} //If an atom is already there, nothing happens
	else{
		//cout<<"addAtom: atom deposited at X: "<<pos[0]<<" and Y: "<<pos[1]<<endl;
		atomList.push_back(Atom (pos,isSticking)); //Copy in tab, if not will be destroyed at the end on function
		updateAtomPointer(atomList.back(),{-1,-1},pos); //-1 because new atom, we set the pointer on the last element of atomList we just push_back()
		return true;
	}


}

double System::energyToProba(double energyBarrier, double dt){
	
	return dt*constPropExp*exp(-energyBarrier/(k_b*T));
}

void System::updateStick(Atom &At){ //set isSticking to true if the atom is beside of that is sticking
	//cout<<"STICK?"<<endl;
	array<int,2> pos=At.get_r();
	//cout<<endl<<endl<<"updateStick: "<<endl;
	//cout<<"--Pos updateStick is X:"<<pos[0]<<" Y: "<<pos[1]<<endl;
	/*
	for(int k(0);k<sizeX;k++){
		for(int l(0);l<sizeY;l++){
			if(fcc[k][l]!=NULL){cout<<"FCC state before defining pointer: X: "<<k<<" Y: "<<l<<" posX: "<<fcc[k][l]->get_r()[0]<<" posY: "<<fcc[k][l]->get_r()[1]<<endl;}
		}
	}*/
	array<bool,6> posNb=posNeighbours(pos);
	/*cout<<"Bool: ";
	for (int m(0);m<6;m++){cout<<posNb[m]<<" ";}
	cout<<endl;*/

	array<array<int,2>,6> indexNb=indexNeighbours(pos);
	//for(int i(0);i<6;i++){cout<<"index "<<i<< " X:"<<indexNb[i][0]<<" Y: "<<indexNb[i][1]<<endl;}
	for (int j(0);j<6;j++){//passe tjr par lapour changer stick
		if(posNb[j]==true){ //if there is a neighbour
			Atom* atomPointer = fcc[indexNb[j][0]][indexNb[j][1]];
			//cout<<"atom pointer "<<j<< " position X: "<< atomPointer ->get_r()[0]<<" and Y: "<<atomPointer ->get_r()[1]<<"and isSticking value: "<<atomPointer -> getSticking()<<endl;
			//if(posNb[j]==true && fcc[indexNb[j][0]][indexNb[j][1]]-> getSticking()){cout<<"STUUUUCK"<<endl;}//Attention au ==tue, important je ne sais pas pourquoi
			if( (atomPointer-> getSticking()==true)){ //PROBLEME ICI: RENTRE PAS Devrait stic les statics avec le update en dehors de isMoving
				At.setSticking(true); //if we are beside a sticky atom, we become sticky too
				//cout<<"---------STICK"<<endl;
				return;
			}
		}
	}

}

array<double,6> System::computeProba(Atom &A){
	array<double,6> result={};

	if(A.getSticking()==false){
		array<bool,6> posNb=posNeighbours(A.get_r());//cout<<"not a stick"<<endl;
		for (int i(0);i<6;i++){
			if (posNb[i]==true){result[i]=0.0;} //if there is already an atom (to wich we don't stick), we can't go there
			else{result[i]=energyToProba(E_M,dt);} //else we do a normal step E_M
		}

	}

	else{//cout<<"atom can't move (stick)"<<endl;
		//return{0.1,0.1,0.1,0.1,0.1,0.1};
		//return{0.0,0.0,0.0,0.0,0.0,0.0};
	
		array<int,2> pos(A.get_r());
		array<array<int,2>,6> index=indexNeighbours(pos);
		int nTouching (numberOfNeighbours(pos)); //How many atoms we are touching
		if(nTouching>=3){return {-1.0,0.0,0.0,0.0,0.0,0.0};}
		else if(nTouching==0){ //This atom is moving freely for now
			for(int i(0);i<6;i++){
				int nTouchingIfMove(numberOfNeighbours(index[i])); //We look if we will be in contact when we move here
				if (nTouchingIfMove==1){result[i]=energyToProba(E_M,dt);} //Nothing around, normal step with energy barrier E_M
				else{result[i]=energyToProba(E_M,dt);} //We suppose the attraction lowers the barrier of a factor 2 if we put /2
			}
		}
		else{ //This atom is already part of a cluster
			array<bool,6> posNb=posNeighbours(pos);
			array<array<int,2>,6> indexNb=indexNeighbours(pos);
			result={0.0,0.0,0.0,0.0,0.0,0.0};//initialize to zero, will fill later
			switch(nTouching){
				case 1: //we are touching a corner
				for (int i(0);i<6;i++){
					if(posNb[i]==true){
						int indexParity=i%2; //if the atom is in even spot (0mod2): +1 is B-step, -1 is A-step. Opposite if odd index
						int spotPlus1= (i+1)%6; //+1 sens aiguilles montre, reste a savoir si A ou B step
						int numberNbSpotPlus1=numberOfNeighbours(indexNb[spotPlus1]);//To know if it is corner-to-corner or corner-to-edge or else
						int spotMinus1= (i+5)%6;//coordinate -1 to avoid negative modulo
						int numberNbSpotMinus1=numberOfNeighbours(indexNb[spotMinus1]);
											

						if (numberNbSpotPlus1==2 && indexParity==0){result[spotPlus1]=energyToProba(E_Bc,dt);} //corner-to-corner B-step for +1
						if (numberNbSpotMinus1==2 && indexParity==0){result[spotMinus1]=energyToProba(E_Ac,dt);} //corner-to-corner A-step for -1
						if (numberNbSpotPlus1==2 && indexParity==1){result[spotPlus1]=energyToProba(E_Ac,dt);} //corner-to-corner A-step for +1
						if (numberNbSpotMinus1==2 && indexParity==1){result[spotMinus1]=energyToProba(E_Bc,dt);} //corner-to-corner B-step for -1
						//Corner to edge: We use edge to corner for now. How to handle special cases >3? for now like corner to edge
						if (numberNbSpotPlus1>2 && indexParity==0){result[spotPlus1]=energyToProba(E_Bc,dt);} //corner-to-edge B-step for +1 	TRY WITH E_AC INSTEAD OF E_AEC
						if (numberNbSpotMinus1>2 && indexParity==0){result[spotMinus1]=energyToProba(E_Ac,dt);} //corner-to-edge A-step for -1
						if (numberNbSpotPlus1>2 && indexParity==1){result[spotPlus1]=energyToProba(E_Ac,dt);} //corner-to-edge A-step for +1
						if (numberNbSpotMinus1>2 && indexParity==1){result[spotMinus1]=energyToProba(E_Bc,dt);} //corner-to-edge B-step for -1


	//OTHER CASES? DETACH?EDGE TO CORNER? FILL A GAP WITH >2 TOUCHING?
					}
				}
				break;
				case 2: //we are touching an edge
				for (int i(0);i<6;i++){
					if (posNb[i]==true && posNb[i+1]==true){//if the two are alligned, else this atom is part of a thin arm we don't want to break
						int indexParity=i%2; //Same as before: (0mod2)-> +1 is B-step, -1 is A-step. Opposite if odd index
						int spotPlus2= (i+2)%6; //since i is the first of the two clockwise, the next free is at +2
						int numberNbSpotPlus2=numberOfNeighbours(indexNb[spotPlus2]);//To know if it is corner-to-corner or corner-to-edge or else
						int spotMinus1= (i+5)%6;//spot -1
						int numberNbSpotMinus1=numberOfNeighbours(indexNb[spotMinus1]);
											

						if (numberNbSpotPlus2==2 && indexParity==0){result[spotPlus2]=energyToProba(E_Aec,dt);} //edge-to-corner: both A-step if 0mod2
						if (numberNbSpotMinus1==2 && indexParity==0){result[spotMinus1]=energyToProba(E_Aec,dt);} //edge-to-corner A-step for -1
						if (numberNbSpotPlus2==2 && indexParity==1){result[spotPlus2]=energyToProba(E_Bec,dt);} //edge-to-corner: both B-step if 1mod2
						if (numberNbSpotMinus1==2 && indexParity==1){result[spotMinus1]=energyToProba(E_Bec,dt);} //edge-to-corner B-step for -1
						//Corner to edge: We use edge to corner for now. How to handle special cases >3? for now like corner to edge
						if (numberNbSpotPlus2>2 && indexParity==0){result[spotPlus2]=energyToProba(E_Ae,dt);} //edge A-step for +2 (for both +2 and -1 if 0mod2)
						if (numberNbSpotMinus1>2 && indexParity==0){result[spotMinus1]=energyToProba(E_Ae,dt);} //edge A-step for -1
						if (numberNbSpotPlus2>2 && indexParity==1){result[spotPlus2]=energyToProba(E_Be,dt);} //edge B-step for +2 (for both +2 and -1 if 1mod2)
						if (numberNbSpotMinus1>2 && indexParity==1){result[spotMinus1]=energyToProba(E_Be,dt);} //edge B-step for -1



					}

				}

				break;
				/*case 3: //can't move from now
				break;
				case 4: //can't move, will return the null array initialized befor the switch
				break;
				case 5: //can't move 
				break;
				case 6: //can't move
				break;*/ //separate if at the begining for optimization
				default:
				cout<<"Error with probability computation, atom won't move"<<endl;
				result={0.0,0.0,0.0,0.0,0.0,0.0};
				break;

			}
		} 
		
	}	
	return result;
	//int num(numberOfNeighbours(pos)); //hit and stick method to test the code
	//if (num!=0){return{0.0,0.0,0.0,0.0,0.0,0.0};} //Bouge pas si se touchent
	//else{return{0.05,0.05,0.05,0.05,0.05,0.05};} //A changer
}

void System::moveSystem(){
	for(size_t i(0);i<atomList.size();i++){ //For all atoms in te list
		
		array<double,6> probaMove (computeProba(atomList[i])); //proba to move at surounding 6 positions

		if (probaMove[0]!=-1.0){ //Optimize simulation speed, if can't move we just return -1 as first probability that will skip the next useless steps
			array<int,2> oldPos(atomList[i].get_r());//get position of the one we are interested
			array<array<int,2>,6> indexNeighb (indexNeighbours(oldPos));
			double randomNumber (distribMove(generator));
			atomList[i].move(probaMove,indexNeighb,randomNumber);
			
			if(atomList[i].isMovingStatus()==true){
				updateAtomPointer(atomList[i],oldPos,atomList[i].get_r()); //if it has moved we update his position in the pointer map
				updateStick(atomList[i]);
			} //and we check if is is sticking to a neighbour atom to update stick status
			//if(atomList[i].getSticking()){cout<<i<<" OK"<<endl;}
			//else{cout<<i<<" NON"<<endl;}
		}
	}
/*
 //We concluded that fcc doesn't contains corrupted pointers
	for(int k(0);k<sizeX;k++){
		for(int l(0);l<sizeY;l++){
			if(fcc[k][l]!=NULL){cout<<"CHECK FCC: X: "<<k<<" Y: "<<l<<" posX: "<<fcc[k][l]->get_r()[0]<<" posY: "<<fcc[k][l]->get_r()[1]<<endl;}
		}
	}

*/

}