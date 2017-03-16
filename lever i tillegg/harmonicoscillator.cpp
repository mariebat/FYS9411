#include "harmonicoscillator.h"
#include <math.h>

#include "hermite.h"

#include <iostream>
using namespace std;


harmonicOscillator::harmonicOscillator()
{
	//Do nothing, default constructur 
}


harmonicOscillator::harmonicOscillator(int index)
{
	//Constructor accepting the index, which maps the variable index to 
	//the quantum numbers (nx, ny, sigma), storing them
	map(index);  //void function, saves nx, ny and sigma (private variables) 
}

harmonicOscillator::harmonicOscillator(int Nx, int Ny, int Sigma)
{
	//Constructor which accepts three quantum numbers, and stores them as private variables. 
	
	nx = Nx; 
	ny = Ny; 
	sigma = Sigma; 

}

int harmonicOscillator::map2p(int Nx, int Ny, int Sigma)
{
	//Accepts the quantum numbers nx, ny and sigma, and computes and outputs the orbital index p 
	
	int p; 

	if (Nx > 9 || Nx < 0)
	{
		cout << "In harmonicOscillator::map2p: Wrong value for quantum number nx entered in map2p. I only accept values from 0 to 9. Aborting program and returning -1.";
		cout << endl;
        return -1;
	}

	if (Ny > 9 || Ny < 0)
	{
		cout << "In harmonicOscillator::map2p: Wrong value for quantum number ny entered in map2p. I only accept values from 0 to 9. Aborting program and returning -1.";
		cout << endl;
        return -1;
	}

	if (Sigma == 1)
	{
		p = (Nx + Ny + 1)*(Nx + Ny) + (2*Ny);
	}
	else if (Sigma == -1)
	{
		p = (Nx + Ny + 1)*(Nx + Ny) + (2*Ny) + 1;
	}
	else
	{
		cout << "In harmonicOscillator::map2p: Wrong value for spin quantum number sigma entered in map2p. I only accept +1 and -1. Aborting program and returning -1 (for p).";
		cout << endl;
        return -1;
	}

	return p;
}

void harmonicOscillator::map(int index)
{
	//Function for mapping orbital indices to quantum numbers nx, ny and sigma
	
	
	//Accepts orbital index p, and calculates the quantum numbers 
	//nx, ny and sigma, storing them within the class private variables. 
	
	int shell; 
	int p = index; 

	if (p == 0)
	{
		nx = 0;
		ny = 0;
		shell = 1; 
		sigma = 1;
        return;
	}
	else if (p == 1)
	{
		nx = 0;
		ny = 0;
		shell = 1; 
		sigma = -1;
        return;
	}
	else if (p > 1 && p < 6)
	{
        shell = 2;
	}
	else if (p > 5 && p < 12)
	{
        shell = 3;
	}
	else if (p > 11 && p < 20)
	{
        shell = 4;
	}
	else if (p > 19 && p < 30)
	{
        shell = 5;
	}
	else if (p > 29 && p < 42)
	{
        shell = 6;
	}
	else if (p > 41 && p < 56)
	{
        shell = 7;
	}
	else if (p > 55 && p < 72)
	{
        shell = 8;
	}
	else if (p > 71 && p < 90)
	{
        shell = 9;
	}
	else if (p > 89 && p < 110)
	{
        shell = 10;
	}
	else //(p < 0 || p > 109)
	{
		cout << "In harmonicOscillator::map: Wrong value for index p entered in map. It must be 0 or larger, and smaller than 110. Aborting program.";
		cout << endl;
		return;
	}

	if (p % 2 == 0)
	{
		sigma = 1;
		nx = (shell - 1) - (p - shell*(shell - 1)) / 2;
		ny = (p - shell*(shell - 1)) / 2;
	}
	else
	{
		sigma = -1;
		nx = (shell - 1) - ((p - 1) - shell*(shell - 1)) / 2;
		ny = ((p - 1) - shell*(shell - 1)) / 2;
	}
}


int harmonicOscillator::map2sigma(int index)
{
	//This function accepts an orbital index between 0 and 109 (shells between 1 and 10), 
	//and computes, stores and outputs the spin quantum number. 
	
	//cannot give shell, nx or ny, only sigma 
	int p = index; 

	if (p < 0 || p > 109)
	{
		cout << "In harmonicOscillator::map2sigma: Wrong index entered in map2sigma. Must be > 0 and < 110. Aborting program and returning 0 for sigma."; 
		cout << endl; 
		return 0; 
	}

	if (p % 2 == 0)
		sigma = 1;
	else
		sigma = -1;

	return sigma;
}

int harmonicOscillator::map2nx(int index)
{ 
	//This function accepts orbital index p and calculates the quantum number nx, 
	//storing it within the class. 
	
	//The function basically only calls void map. 
	
	//If the given index is out of bounds (only shells 1 to 10 accepted), 
	//map outputs an error message and exits. 
	
	map(index); 

	return nx; 
}

int harmonicOscillator::map2ny(int index)
{
	//This function accepts orbital index p and calculates the quantum number ny, 
	//storing it within the class. 
	
	//The function basically only calls void map. 
	
	//If the given index is out of bounds (only shells 1 to 10 accepted), 
	//map outputs an error message and exits. 

	map(index); 

	return ny; 
}

double harmonicOscillator::wavefunction(int index, double x, double y, double omega)
{
	//This function computes the Harmonic oscillator wavefunction phi at spatial coordinates x and y,
	//for specific quantum numbers nx and ny. The spin quantum number is not relevant. The harmonic oscillator frequency must be specified. 

    double phi;
    double Hx;
    double Hy;
	
	//Find quantum numbers nx and ny: 
	map(index); 

	//Copmpute amplitudes
    double Ax = pow((omega/pi),0.25)*(1/sqrt(pow(2,nx)*factorial(nx)));
    double Ay = pow((omega/pi),0.25)*(1/sqrt(pow(2,ny)*factorial(ny)));

	//Compute polynomials 
	Hx = hermPol.calcHermite(nx, sqrt(omega)*x); 
	Hy = hermPol.calcHermite(ny, sqrt(omega)*y); 

	
	//Compute wavefunction  
    phi = Ax*Ay*Hx*Hy*exp(-omega*0.5*(x*x + y*y));

	return phi; 
}

int harmonicOscillator::factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}



int harmonicOscillator::energy(int index)
{
	//Harmonic oscillator energy as a function of orbital index 
	int energy; 
	map(index); 

	energy = nx + ny + 1; 
	return energy; 
}

int harmonicOscillator::magicNumber(int index)
{
	int N; 
	map(index); 

	N = (nx + ny + 1)*(nx + ny + 2);
	return N; 
}

int harmonicOscillator::degeneracy(int index)
{
	int deg; 
	map(index); 

	deg = 2*(nx + ny + 1 ); 

	return deg; 
}

void harmonicOscillator::print()
{
	//This function prints the quantum numbers (nx, ny, sigma) which are currently 
	//stored as private variables. 

    cout << "nx is: ";
    cout << nx;
    cout << ".\n";

    cout << "ny is: ";
    cout << ny;
    cout << ".\n";

    cout << "sigma is: ";
    cout << sigma;
    cout << ".\n";
}
