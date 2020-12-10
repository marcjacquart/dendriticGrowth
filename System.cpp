#include "System.h"
#include<cmath>
using namespace std;

//Methods for the System Class:

void System::updateAtomPointer(Atom &A, array<int,2>const (&oldPos), array<int,2> const (&newPos)){ // Update the fcc pointer 2d array
	
	if(oldPos[0]!=-1){fcc[oldPos[0]][oldPos[1]]=NULL;}		// Destroy old pointer, -1 means that the atom is newly created, no need to destroy its old position.
	
	fcc[newPos[0]][newPos[1]]=&A; 							// Assign in tab the adress of A, wich is a reference (just another name for the variable) over the atom we want
	
}


void System::printSystem(string const &outFile, bool printAdatom){								// Print the simulation results: position of the atoms in the system
	ofstream writeTxt;																			// Define output file stream
	writeTxt.open(outFile,ios_base::app); 														// Text documment with results
	if (not writeTxt.fail()){
		for(int i(0);i<nAtomsInSyst;i++){
			if(atomList[i].getSticking()==true||printAdatom==true){
				double coordX(atomList[i].get_r()[0]);
				double coordY(atomList[i].get_r()[1]);
				writeTxt<< fmod(coordX+0.5*coordY,sizeX)<<","<< fmod(0.866*coordY,sizeY)<<";"; 	// Gives: x_0,y_0;x_1,y_1;x_2,y_2;... fmod:modulo for double
			}
		}
		writeTxt<<endl; 																		// One line for each timestep
	}
	writeTxt.close();
}

void System::printParam(string const &paramFile){												// Print physical parameters in the text file at the begining
	ofstream writeParam;
	writeParam.open(paramFile,ios_base::app);
	writeParam<<"dt= "<<get_dt()<<endl;
	writeParam<<"Number of steps: "<<100/dt<<endl;
	array<string,7> energyTab={"E_M","E_Ac","E_Bc","E_Ae","E_Be","E_Aec","E_Bec"};
	for(int j(0);j<7;j++){writeParam<<"Probability movement "<<energyTab[j]<<" = "<<energyToProbaTab[j]<<endl;}

	writeParam.close();

}


array<array<int,2>,6> System::indexNeighbours(array<int,2> const (&pos)){						// index {X,Y} of each neighbour of one given Atom
	
	array<int,2>index1={pos[0],(pos[1]+1)%sizeY};
	array<int,2>index2={(pos[0]+1)% sizeX,pos[1]};
	array<int,2>index3={(pos[0]+1)%sizeX,(pos[1]-1+sizeY)%sizeY}; 								// +sizeY to avoid -1%sizeY=-1 :/ reminder not modulo (problem for negative values)
	array<int,2>index4={pos[0],(pos[1]-1+sizeY)%sizeY};
	array<int,2>index5={(pos[0]-1+sizeX)%sizeX,pos[1]};
	array<int,2>index6={(pos[0]-1+sizeX)%sizeX,(pos[1]+1)%sizeY};
	return{index1,index2,index3,index4,index5,index6};

}


int System::numberOfNeighbours (array<int,2>const (&pos)){ 		// Check if there is something in each of the 6 surounding spots to return the  number of neighbouring atoms of one site
	array<array<int,2>,6> index =indexNeighbours(pos);
	int result(0);
	for(int i(0);i<6;i++){
		if(fcc[index[i][0]][index[i][1]]!=NULL){result+=1;}
	}
	return result;

}
array<bool,6> System::posNeighbours(array<int,2>const (&pos)){ // Index 1 to 6 with respect to convention (see README.txt), true if an atom is here
	array<bool,6> results ={false,false,false,false,false,false};
	array<array<int,2>,6> index =indexNeighbours(pos);
	for(int i(0);i<6;i++){
		if(fcc[index[i][0]][index[i][1]]!=NULL){results[i]=true;}
	}
	return results;
}

bool System::addAtom(array<int,2> const (&pos), bool isSticking){ 				// Add an atom to the system, not sticking by default
	if (fcc[pos[0]][pos[1]]!=NULL || nAtomsInSyst>=maxAtoms){return false;} 	// If an atom is already there or we are at max number of allowed atoms, nothing happens
	else{
		Atom A(pos,isSticking);													// Create the new atom
		atomList[nAtomsInSyst]=A; 												// Copy in tab, if not will be destroyed at the end on function. 
																				// /!\ Do not use push back, rearrangement in memory corrupts the fcc pointer array
		updateAtomPointer(atomList[nAtomsInSyst],{-1,-1},pos); 					// -1 because new atom (for speed optimization in updateAtomPointer), 
																				// we set the pointer on the last element of atomList we just added to the atomList
		nAtomsInSyst=nAtomsInSyst+1;											// Update the number of atoms in the system
		
		return true;
	}


}


void System::updateStick(Atom &At){ 											// Set isSticking to true if the atom is beside of that is sticking
	array<int,2> pos=At.get_r();

	array<bool,6> posNb=posNeighbours(pos);

	array<array<int,2>,6> indexNb=indexNeighbours(pos);
	for (int j(0);j<6;j++){
		if(posNb[j]==true){ 													// If there is a neighbour
			Atom* atomPointer = fcc[indexNb[j][0]][indexNb[j][1]];
			if( (atomPointer-> getSticking()==true)){							// Check if it is sticking
				At.setSticking(true); 											// If we are beside a sticky atom, we become sticky too
				return;															// Once we setSticking(true), no point to continue the loop
			}
		}
	}

}

double System::energyToProba(double energyBarrier, double dt){ 					// Used at the begining to compute the probability of movement from the energy barriers
	
	return dt*constPropExp*exp(-energyBarrier/(k_b*T));
}

array<double,6> System::computeProba(Atom const &A){							// Compute the probability movement of an atom for the 6 surounding positions
	array<double,6> result={0.0,0.0,0.0,0.0,0.0,0.0};							// Will add the results =/= 0 afterwards

	if(A.getSticking()==false){ 												// This part is used only if we chose the cluster position by generatimg non sticking atoms
		array<bool,6> posNb=posNeighbours(A.get_r());
		for (int i(0);i<6;i++){
			if (posNb[i]==false){result[i]=energyToProbaTab[0];} 				// If nobody there we do a normal step, else if an atom is already there the probability stay 0.0 as initialized
			
		}
		
		return result;
	}

	else{
		//return{0.0,0.0,0.0,0.0,0.0,0.0}; 										// /!\ Uncomment this line for hit and stick: wont go further down, if stick just dont move at all, work for one big cluster (generated particles don't stick at first)
		array<int,2> pos(A.get_r());
		
		int nTouching (numberOfNeighbours(pos)); 								// How many atoms we are touching
		if(nTouching>=2){return {-1.0,0.0,0.0,0.0,0.0,0.0};}					// Neglect >1 to gain speed at low energy, -1.0 for time optimization in moveSystem()
																				// (flag to tell other function that this atom doesn't move, no need to compute anything or draw a random number for this one)
		
		else if(nTouching==0){ 													// This atom is moving freely for now (touches nobody)

			double resFree=energyToProbaTab[0];									// Probability tab for not computing each time with energyToProba, we only have 7 potentially different energy barriers, 
																				// compute all at the begining and the just acces the tab where the result are stored.
			return{resFree,resFree,resFree,resFree,resFree,resFree};
			/* 																	// Easy approach effect, not treated here. Uncomment to activate
			array<array<int,2>,6> index=indexNeighbours(pos);
			for(int i(0);i<6;i++){
				int nTouchingIfMove(numberOfNeighbours(index[i])); 				// We look if we will be in contact when we move here
				if (nTouchingIfMove==1){result[i]=energyToProba(E_M,dt);} 		// Nothing around, normal step with energy barrier E_M
				else{result[i]=energyToProba(E_M,dt);} 							// We suppose the attraction lowers the barrier of a factor 2 if we put /2
			}*/ 
		}
		else{																	// This atom is already part of a cluster, with simplified move this means that the atom touches only one other atom
			array<bool,6> posNb=posNeighbours(pos);
			array<array<int,2>,6> indexNb=indexNeighbours(pos);
			switch(nTouching){
				case 1: 														// We are touching a corner
				for (int i(0);i<6;i++){
					if(posNb[i]==true){
						int indexParity=i%2; 									// If the atom is in even spot (0mod2): +1 is B-step, -1 is A-step. Opposite if odd index
						int spotPlus1= (i+1)%6; 								// +1 clockwise, now just need to know if it is an A or B-step
						int numberNbSpotPlus1=numberOfNeighbours(indexNb[spotPlus1]); 	// To know if it is corner-to-corner or corner-to-edge or else
						int spotMinus1= (i+5)%6;										// Coordinate -1 to avoid negative modulo
						int numberNbSpotMinus1=numberOfNeighbours(indexNb[spotMinus1]);
											
						// Faster with the tab (calculated once for all), E_Bc neglected
						// Uncomment to allows movements
						//if (numberNbSpotPlus1==2 && indexParity==0){result[spotPlus1]=energyToProba(E_Bc,dt);} 	// Corner-to-corner B-step for +1
						if (numberNbSpotMinus1==2 && indexParity==0){result[spotMinus1]=energyToProbaTab[1];}		// Corner-to-corner A-step for -1
						if (numberNbSpotPlus1==2 && indexParity==1){result[spotPlus1]=energyToProbaTab[1];} 		// Corner-to-corner A-step for +1
						//if (numberNbSpotMinus1==2 && indexParity==1){result[spotMinus1]=energyToProba(E_Bc,dt);} 	// Corner-to-corner B-step for -1

						//Corner to edge: We use edge to corner for now. How to handle special cases >3? for now like corner to edge
						//if (numberNbSpotPlus1>2 && indexParity==0){result[spotPlus1]=energyToProba(E_Bc,dt);} 	// Corner-to-edge B-step for +1 
						if (numberNbSpotMinus1>2 && indexParity==0){result[spotMinus1]=energyToProbaTab[1];} 		// Corner-to-edge A-step for -1
						if (numberNbSpotPlus1>2 && indexParity==1){result[spotPlus1]=energyToProbaTab[1];} 			// Corner-to-edge A-step for +1
						//if (numberNbSpotMinus1>2 && indexParity==1){result[spotMinus1]=energyToProba(E_Bc,dt);} 	// Corner-to-edge B-step for -1



					}
				}
				break;
				/* 																									// Edge diffusion are disabled, uncomment to enable
				case 2: //we are touching an edge
				for (int i(0);i<6;i++){
					if (posNb[i]==true && posNb[i+1]==true){														// If the two are alligned, else this atom is part of a thin arm we don't want to break
						int indexParity=i%2; 																		// Same as before: (0mod2)-> +1 is B-step, -1 is A-step. Opposite if odd index
						int spotPlus2= (i+2)%6; 																	// Since i is the first of the two clockwise, the next free is at +2
						int numberNbSpotPlus2=numberOfNeighbours(indexNb[spotPlus2]);								// To know if it is corner-to-corner or corner-to-edge or else
						int spotMinus1= (i+5)%6;//spot -1
						int numberNbSpotMinus1=numberOfNeighbours(indexNb[spotMinus1]);
											

						if (numberNbSpotPlus2==2 && indexParity==0){result[spotPlus2]=energyToProba(E_Aec,dt);} 	// Edge-to-corner: both A-step if 0mod2
						if (numberNbSpotMinus1==2 && indexParity==0){result[spotMinus1]=energyToProba(E_Aec,dt);} 	// Edge-to-corner A-step for -1
						if (numberNbSpotPlus2==2 && indexParity==1){result[spotPlus2]=energyToProba(E_Bec,dt);} 	// Edge-to-corner: both B-step if 1mod2
						if (numberNbSpotMinus1==2 && indexParity==1){result[spotMinus1]=energyToProba(E_Bec,dt);} 	// Edge-to-corner B-step for -1

						//Corner to edge: We use edge to corner for now. How to handle special cases >3? for now like corner to edge
						if (numberNbSpotPlus2>2 && indexParity==0){result[spotPlus2]=energyToProba(E_Ae,dt);} 		// Edge A-step for +2 (for both +2 and -1 if 0mod2)
						if (numberNbSpotMinus1>2 && indexParity==0){result[spotMinus1]=energyToProba(E_Ae,dt);} 	// Edge A-step for -1
						if (numberNbSpotPlus2>2 && indexParity==1){result[spotPlus2]=energyToProba(E_Be,dt);} 		// Edge B-step for +2 (for both +2 and -1 if 1mod2)
						if (numberNbSpotMinus1>2 && indexParity==1){result[spotMinus1]=energyToProba(E_Be,dt);} 	// Edge B-step for -1



					}

				}

				break;*/ 		//Desactivated for speed, negligable
				/*case 3: 		// Can't move from now
				break;
				case 4: 		// Can't move, will return the null array initialized befor the switch
				break;
				case 5: 		// Can't move 
				break;
				case 6: 		// Can't physically move
				break;*/ 		
				default:
				cout<<"Error with probability computation, atom won't move"<<endl; // Will return the 0.0 tab, print an error message
				break;

			}
		} 
		
	}	
	return result;				// Return the filled probability tab
}


void System::printFccStatus(){	// For debug: print all the position of all atoms in the fcc pointer tab (to verify if they are not corrupted)
	for(int i(0);i<sizeX;i++){
		for(int j(0);j<sizeY;j++){
			if(fcc[i][j]!=NULL){cout<<"X: "<<fcc[i][j]->get_r()[0]<<" Y: "<<fcc[i][j]->get_r()[1]<<", i: "<<i<<" j: "<<j<<endl;}
		}
	}
	cout<<endl;
}

void System::moveSystem(){													// Move every atom in the atomList

	for(int i(0);i<nAtomsInSyst;i++){ 										// For every atoms in the list
		if(atomList[i].getSticking()==false){updateStick(atomList[i]);} 	// If not sticking and not connected to a cluster yet, we verify that it is not sticking before moving
		array<double,6> probaMove (computeProba(atomList[i])); 				// Probability to move at surounding 6 positions

		if (probaMove[0]!=-1.0){ 											// Optimize simulation speed, if can't move we just return -1 as first probability that will skip the next useless steps
			array<int,2> oldPos(atomList[i].get_r());						// Get position of the one we are interested
			array<array<int,2>,6> indexNeighb (indexNeighbours(oldPos));
			double randomNumber (distribMove(generator));					// Generate a random number to select the movement
			atomList[i].move(probaMove,indexNeighb,randomNumber);
			
			if(atomList[i].isMovingStatus()==true){
				updateAtomPointer(atomList[i],oldPos,atomList[i].get_r()); 	// If it has moved we update his position in the pointer map
				
			} 
		}
	}
}