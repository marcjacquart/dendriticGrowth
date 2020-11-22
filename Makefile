#Makefile for the Simulation Project: From tutorial https://www.cs.bu.edu/teaching/cpp/writing-makefiles/

#Variables:
CXX =g++
CXXFLAGS = -std=c++11 -Wall -g 


#Targets:
main: main.o System.o Atom.o
	$(CXX) $(CXXFLAGS) -o main main.o System.o Atom.o
	
	
main.o: main.cpp System.h Atom.h
	$(CXX) $(CXXFLAGS) -c main.cpp
	
System.o: System.cpp System.h Atom.h
	$(CXX) $(CXXFLAGS) -c System.cpp

Atom.o: Atom.cpp Atom.h
	$(CXX) $(CXXFLAGS) -c Atom.cpp
	
#Cmd qu'on aimerait que ca fasse:	
# g++ -c main.cpp
# g++ -c Atom.cpp
# g++ -c System.cpp

#g++ -o main main.o Atom.o System.o
