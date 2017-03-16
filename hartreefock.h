#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include "coulombint.h"
#include <armadillo>

using namespace arma;

//The Hartree-Fock class is defined here. 

//The class contains all the functions necessary for performing a self-consistent Hartree-Fock calculation for a 
// 2D all-electron quantum dot system. The electrons reside in a harmonic oscillator potential, and a harmonic 
//oscillator basis is used as the starting point for the wavefunction. 

//The claculations are performed in polar coordinates. 

//The C++ library Armadillo is employed for matrix, vector and eigenvalue operations. 

//The two-body interactions are computed using the program supplied by Morten Hjorth-Jensen. 
//These functions are stored within the class coulombInt. 

class HartreeFock
{
public:

    int factorial(int n); //Just in case it's needed

    //Index and quantum number mapping functions
    void map(int p);

    int map2n(int index); //calls map
    int map2m (int index); // calls map
    int map2sigma(int index); //Basically only checks whether p is even or odd 

    int map2p(int a, int b, int s); //finds orbital index depending on quantum numbers n, m and sigma 


    //Functions accepting quantum numbers and calculating harmonic oscillator properties 
    int HOenergy(int index); 
    int magicNumber(int index);
    int degeneracy(int index);
    
	
	//Final part: Hartree-Fock code! 
	
	//Calculating the Coloumb interaction: 
    double computeCoulombElement(double hw, int p, int q, int r, int s); //Computes one Coulomb two-body matrix element (on the fly)
	//vec getTwoBodyVector(double hw, int nOrbitals); //Not made; should store two-body elements in a vector 
	
	//Functions to get necessary matrices: 
    mat getDensityMatrix(mat C, int nOrbitals, int nElectrons); //Computes and outputs the density matrix rho (which is the product of columns in C)
	mat compHOenergyMatrix(double hw, int nOrbitals); //only need this if matrix elements are prestored 
    mat compFockMatrix(mat C, int nOrbitals, int nElectrons, double hw); //makes the Fock matrix, the hamiltonian h_HF, which is diagonalized in doSCL()
    
	//Function to compute the total Hartree-Fock energy:
    double compHFenergy(mat C, vec HFenergies, int nOrbitals, int nElectrons, double hw); 

	//Function to do the self-consistent loop:
    void doSCL(double hw, int maxLoops, int nOrbitals, int nElectrons); 


private:
    // Coulomb interaction (relevant functions are stored in this class)
    coulombInt Coulomb;

    //quantum numbers
    int n;
    int m;
    int sigma;


};

#endif // HARTREEFOCK_H
