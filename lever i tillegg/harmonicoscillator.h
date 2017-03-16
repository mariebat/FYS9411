#include "hermite.h"


#ifndef HARMONICOSCILLATOR_H
#define HARMONICOSCILLATOR_H

//This class defines the harmonic oscillator basis to be used in Hartree-Fock calculations. 

//The basis is constructed in Cartesian coordinates, and was therefore not used in the 
//Hartree-Fock calculations for Project 1. Parts of the class might be used in Project 2. 

//The main capabilities of the class are to map between orbital indices and quantum numbers, 
//and to compute the wavefunction for a specific orbital, oscillator frequency and spatial 
//position. 

//The harmonic oscillator class also uses the hermite polynomial class. 

class harmonicOscillator
{
public:
	harmonicOscillator(); 
	harmonicOscillator(int index); 
	harmonicOscillator(int Nx, int Ny, int Sigma); 

	int map2p(int Nx, int Ny, int Sigma); 

	void map(int index); 

	int map2sigma(int index);
	int map2nx(int index); 
	int map2ny(int index); 
	
    double wavefunction(int index, double x, double y, double omega);

	int energy(int index); 
	int magicNumber(int index);
	int degeneracy(int index); 

    int factorial(int n);

    void print();

private:
	int nx;
	int ny;
	int sigma;

    double pi = 3.1415;

	hermite hermPol;

};

#endif // HARMONICOSCILLATOR_H

