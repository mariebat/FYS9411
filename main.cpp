#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <time.h>

#include "hartreefock.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{

	cout << "----HF solver----" << endl;
	cout << "This program considers a 2D quantum dot system. " << endl;
	cout << "Need to specify the frequency, the number of iterations the program can do, the number of orbitals included (cut-off) and the number of electrons included. " << endl;
    
    //Initialize timer
    clock_t start, stop;


	//Initialize Hartree-Fock class
	HartreeFock HF;
	
	//Define necessary variables 
    double hw = 0.1; //harmonic oscillator frequency
	
    int nShells = 5; //number of shells (cut-off)
    int nOrbitals = nShells*(nShells+1); //the number of single particle states included
    int nElectrons = 12; //the number of electrons included in the calculation
    //The system must be closed-shell! That is, #electrons = fermi level, and nElectrons must be 2, 6, 12, 20, 30, 42, 56, 72, 90, or 110.
	
	//Decide the number of times the program can try to achieve self-consistency 
    int maxLoops = 100;
	
    start = clock();
    //run the self-consistent loop
    HF.doSCL(hw, maxLoops, nOrbitals, nElectrons);
    stop = clock();

    cout << "Time used:" << (stop-start)/CLOCKS_PER_SEC << "s." << endl;

    return 0;
}


