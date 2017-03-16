#include <math.h>
#include <iostream>
#include <armadillo>
#include <iomanip>
#include "hartreefock.h"
using namespace std;
using namespace arma;

//The member functions of the Hartree-Fock class is defined here. 

//The class contains all the functions necessary for performing a self-consistent Hartree-Fock calculation for a 
// 2D all-electron quantum dot system. The electrons reside in a harmonic oscillator potential, and a harmonic 
//oscillator basis is used as the starting point for the wavefunction. 

int HartreeFock::factorial(int a)
{
  return (a == 1 || a == 0) ? 1 : factorial(a - 1) * a;
}


//Functions accepting quantum numbers and calculating harmonic oscillator properties::

int HartreeFock::HOenergy(int index) //Accepts orbital index, and finds harmonic oscillator energies 
{
    map(index);

    int energy = 2*n + abs(m) + 1;
    return energy;
}

int HartreeFock::magicNumber(int index) //Accepts orbital index, and finds the quantum dot magic number 
{
    map(index);

    int N = (2*n + abs(m) + 1)*(2*n + abs(m) + 2);
    return N;
}

int HartreeFock::degeneracy(int index) //Accepts orbital index, and finds degeneracy of given orbital 
{
    map(index); 

    int deg = 2*(2*n + abs(m) + 1 );

    return deg;
}


//Index and quantum number mapping functions:: 

void HartreeFock::map(int p)
{
	//This function accepts orbital index p, and calculates the quantum numbers
    //n, m and sigma, storing them as class private variables.
	
    int shell;


    if (p == 0)
    {
        n = 0;
        m = 0;
        sigma = 1;
        return;
    }
    else if (p == 1)
    {
        n = 0;
        m = 0;
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
        cout << "In HartreeFock::map: Wrong value for index p entered in map. It must be 0 or larger, and smaller than 110. Aborting program.";
        cout << endl;
        return;
    }

	int A = p - shell*(shell - 1); 
	
    if (p % 2 == 0 )
    {
        sigma = 1;
		m = -(shell- 1) + A;  
    }
    else
    {
        sigma = -1;
		m = -shell + A;  
    }

    n = (shell - abs(m) - 1)/2;

}

int HartreeFock::map2n(int index) //Calls map (which checks p not out of bounds) and returns n
{
    map(index);
    return n;
}

int HartreeFock::map2m(int index) //Calls map (which checks p not out of bounds) and returns m
{
    map(index);
    return m;
}

int HartreeFock::map2sigma(int index)  //Checks whether the orbital index is even or odd, and returns the spin
{

	//Could include the check below, but it is redundant; sigma is never computed unless  n and m are computed for the same index at the same time. 
	//These both call map, and here the check for whether p is out of bounds is performed. It is unnecessary to check three times! 

	//if (index < 0 || index > 109)
	//{
	//	cout << "In HartreeFock::map2sigma: Wrong index entered in map2sigma. Must be >= 0 and < 110. Aborting program and returning 0 for sigma."; 
	//	cout << endl; 
	//	return 0; 
	//}
	
	
    if (index % 2 == 0 )
    {
        sigma = 1;
    }
    else
    {
        sigma = -1;
    }

    return sigma;
}

//Don't really need this function, only for testing 
int HartreeFock::map2p(int a, int b, int s) //Accepts n, m and sigma, finds and returns the orbital index p. If QNs are out of bounds, p is returned as -1. 
{
	if (a < 0) 
	{
		cout << "In HartreeFock::map2p(n,m,sigma): Quantum number n cannot be negative. Try again! Aborting program and returning -1 for p." << endl; 
		return -1; 
	} 
    n = a;
    m = b;
    sigma = s;

    int shell = (2*n + abs(m) + 1);
	int A; 
	int p; 

    if (sigma == 1)
    {
        A = m + shell - 1; 
    }
    else if (sigma == -1)
    {
        A = m + shell; 
    }
    else
    {
        cout << "In HartreeFock::map2p: Wrong value for spin quantum number sigma entered in map2p. I only accept +1 and -1. Aborting program and returning -1 for p.";
        cout << endl;
        return -1;
    }

	p = A + shell*(shell-1); 
	
	
	//Shouldn't need this check, if everything else is coded right 
    //if (p < 0)
    //{
    //    cout << "In HartreeFock::map2p: Somehow, p is negative. You messed up. Try again." << endl;
    //    return -1;
    //}
	
    return p;
}

//Final part: Hartree-Fock code! 
	
//Calculating the Coloumb interaction: 
double HartreeFock::computeCoulombElement(double hw, int p, int q, int r, int s) //Input frequency, and 4 orbital indices
{
    double V;
    double V0 = 0.0;
    double V1 = 0.0;

	//Compute all the quantum numbers 
    int np = map2n(p);
    int mp = map2m(p);
    int sp = map2sigma(p);

    int nq = map2n(q);
    int mq = map2m(q);
    int sq = map2sigma(q);

    int nr = map2n(r);
    int mr = map2m(r);
    int sr = map2sigma(r);

    int ns = map2n(s);
    int ms = map2m(s);
    int ss = map2sigma(s);

	
	//Calculate the two-body interaction, direct and exchange, while only connecting the same spins 
    if ((sp == sr) && (sq == ss) )
    {
        V0 = Coulomb.Coulomb_HO(hw, np, mp, nq, mq, nr, mr, ns, ms);
    }

    if ((sp == ss) && (sq == sr))
    {
        V1 = Coulomb.Coulomb_HO(hw, np, mp, nq, mq, ns, ms, nr, mr);
    }

    V = V0 - V1;
    return V;
}

//Functions to get the necessary matrices: 

mat HartreeFock::getDensityMatrix(mat C, int nOrbitals, int nElectrons) //Computes and outputs the density matrix rho (which is the product of columns in C)
{
    mat rho(nOrbitals,nOrbitals,fill::zeros);

	double rsum; 
	
    for (int gamma = 0; gamma < nOrbitals; gamma++)
	{
		for (int delta = 0; delta < nOrbitals; delta++)
		{
			rsum = 0.0; 
			
			for (int i = 0; i < nElectrons; i++)
			{
                rsum += C(i,gamma)*C(i,delta);
			}
			rho(gamma,delta) = rsum; 
		}
		
	}

    return rho;
}

mat HartreeFock::compHOenergyMatrix(double hw, int nOrbitals) //Makes the one-body interaction matrix (harmonic oscillator energies). Only need this if matrix elements are prestored. 
{
	//Only use this function if matrix elements are prestored! 
	
	mat E(nOrbitals,nOrbitals,fill::zeros); 
	
	
	for (int i = 0; i < nOrbitals; i++)
	{
		for (int j = 0; j < nOrbitals; j++)
		{
            if (i == j)
            {
                E(i,j) = hw*HOenergy(i); 
            }

     	}
	}
		
	return E; 
}




mat HartreeFock::compFockMatrix(mat C, int nOrbitals, int nElectrons, double hw) //Makes the Fock matrix, the hamiltonian h_HF, which is diagonalized in doSCL(). 
{
	mat hF(nOrbitals,nOrbitals,fill::zeros); 
	
	//Need this if you prestore matrix elements -->
	//mat Ftemp(nOrbitals,nOrbitals,fill::zeros); 
	//mat E = compEnergyMatrix(int nOrbitals); 
	
	//Get density matrix 
	mat densityMatrix = getDensityMatrix(C, nOrbitals, nElectrons);
	
	double sumhF; 
		
		for (int alpha = 0; alpha < nOrbitals; alpha++)
		{
			for (int beta = 0; beta < nOrbitals; beta++)
			{
                sumhF = 0.0;
                //Two-body interaction (computing matrix elements on the fly)
                for (int gamma = 0; gamma < nOrbitals; gamma++)
                {
                    for (int delta = 0; delta < nOrbitals; delta++)
                    {
                        sumhF += densityMatrix(gamma,delta)*computeCoulombElement(hw, alpha, gamma, beta, delta);
                    }
                }
                hF(alpha,beta) = sumhF;
		
				//One-body interaction; harmonic oscillator energies! Only for diagonal entries
				if (alpha == beta)
				{ 
					hF(alpha,alpha) +=hw*HOenergy(alpha);
				}
				
			}
		}
	
	return hF; 
	
}

//Function to compute the total Hartree-Fock energy:
double HartreeFock::compHFenergy(mat C, vec SPenergies, int nOrbitals, int nElectrons, double hw) //gives you the final energy
{

    //formula v01
    //double Ehf = 0.0;
	double Ehfrho = 0.0;  
	double V;
	
	
	//Get density matrix to speed things upp - a lot! 
	mat densityMatrix = getDensityMatrix(C, nOrbitals, nElectrons);
	
    //formula V02
    //double Ehfalt = 0.0;

    //formula v02: The one-body term:
    //for (int i = 0; i < nElectrons; i++){
    //    for (int alpha = 0; alpha < nOrbitals; alpha++){
    //        Ehfalt += C(i,alpha)*C(i,alpha)*hw*HOenergy(alpha);
    //    }
    //}

    //HF single-particle energies: 
    for (int i = 0; i < nElectrons; i++)
    {
        //Ehf += SPenergies(i);
		Ehfrho += SPenergies(i);
    }

	// The two-body term (with density matrix) - faster! 
	for (int alpha = 0; alpha < nOrbitals; alpha++){
        for (int beta = 0; beta < nOrbitals; beta++){
            for (int gamma = 0; gamma < nOrbitals; gamma++){
                for (int delta = 0; delta < nOrbitals; delta++){
                    V = computeCoulombElement(hw, alpha, gamma, beta, delta);
                    Ehfrho -= 0.5*densityMatrix(alpha,beta)*densityMatrix(gamma,delta)*V;
                }
            }

        }
    }
	

    // The two-body term (without density matrix):
    //for (int i = 0; i < nElectrons; i++){
    //    for(int j = 0; j < nElectrons; j++){
    //
    //        for (int alpha = 0; alpha < nOrbitals; alpha++){
    //            for (int beta = 0; beta < nOrbitals; beta++){
    //                for (int gamma = 0; gamma < nOrbitals; gamma++){
    //                    for (int delta = 0; delta < nOrbitals; delta++){
    //                        V = computeCoulombElement(hw, alpha, gamma, beta, delta);
    //                        //Ehfalt += 0.5*C(i,alpha)*C(j,gamma)*C(i,beta)*C(j,delta)*V;
    //                        Ehf -= 0.5*C(i,alpha)*C(j,gamma)*C(i,beta)*C(j,delta)*V;
    //                    }
    //                }
    //
    //            }
    //        }
    //    }
    //}

    //cout << setprecision(8) << "In compHFenergy; Ehf (with SPenergies): " << Ehf << "." << endl;
	cout << setprecision(8) << "In compHFenergy; Ehf using rho (with SPenergies): " << Ehfrho << "." << endl;
	
    //cout << setprecision(8) << "In compHFenergy; Ehf alternative (with HOenergies):" << Ehfalt << ". " << endl;
    return Ehfrho;
}

//Function to do the self-consistent loop:
void HartreeFock::doSCL(double hw, int maxLoops, int nOrbitals, int nElectrons)
{
    //Check that the system has correct number of orbitals
    if ((nOrbitals < 0 )|| (nOrbitals > 109) )
    {
        cout << "Orbital cut-off is too large or small! Needs to be smaller than 109 (and larger than zero). " << endl;
        cout << "Specify fewer (but positive) shells!" << endl;
        return;
	}
	
	//initialize
    double sum;
	
	//define loop conditions
    int nIterations = 0;
    double treshold = 1.0e-6;
    double difference = 1.0;
	
	//initialize matrices
    mat C(nOrbitals,nOrbitals,fill::eye);
	mat HFmatrix;
	
	
	//Initialize vectors 
	vec SPenergies(nOrbitals,fill::zeros);
	vec oldenergies(nOrbitals,fill::zeros);
	vec newenergies(nOrbitals,fill::zeros);
    //vec diff;
	

    while (nIterations < maxLoops && difference > treshold)
    {
		cout << "Iteration: " << nIterations + 1 << endl; 
		
		//Find Fock matrix for given C (starts with unitary matrix): 
		
		//The function compFockMatrix() initializes mat HFmatrix = zeros(nOrbitals,nOrbitals)
		//It first finds the density matrix, then iterates over four indices, and computes the Fock hamiltonian matrix 
        HFmatrix = compFockMatrix(C,nOrbitals,nElectrons,hw);
		
		//Computing the HF eigenvalue problem
        eig_sym(SPenergies,C,HFmatrix); //updates SPenergies, C, HFmatrix 
        C = trans(C);

		//Check whether single-particle energy vector is consistent with the previous iteration 
		newenergies = SPenergies; 
		
		//Brute force computation of the difference between new and old SP energies 
        sum = 0.0;
        for (int i = 0; i < nOrbitals; i++)
        {
        	sum +=  (abs(newenergies(i)-oldenergies(i)))/nOrbitals;
        }
        difference = sum;
		
		//Alternative computation of difference between new and old SP energies -->
		
        //diff = newenergies - oldenergies;
        //difference = abs(diff.max());
		
        oldenergies = newenergies;
		
        nIterations++;
    }

	//Calculate total HF energy: (time-consuming!)
	double Ehf = compHFenergy(C, SPenergies, nOrbitals, nElectrons,hw); 
	
	//Print results to screen: 
	cout << setprecision(8) << "Final Hartree Fock energy E [a.u.] = " << Ehf <<  endl;
    if (nIterations == maxLoops)
    {
        cout << "Self-consistency not reached. " << endl;
        cout << "Error/ Difference: " << difference << endl;
        for (int i = 0; i < nElectrons; i++)
        {
            cout << "SPenergy nr " << i + 1 << " is " << SPenergies(i) << endl;
        }
    }
    else
    {
        cout << "Self-consistency was reached after " << nIterations << " iterations. " << endl;
        cout << "Allowed error: "  << treshold << endl;
        cout << "Actual difference in SP energies: " << difference << endl;
        cout << "Single particle energies: " << endl;

        for (int i = 0; i < nOrbitals; i++)
        {
            cout << "SPenergy nr " << i + 1 << " is " << SPenergies(i) << endl;
        }

        cout << setprecision(10) << "Final Hartree Fock energy E [a.u.] = " << Ehf <<  endl;
    }
	
	//Look at and save matrices: 
	
    //HFmatrix.print();
    //C.print();
    //C.save("C_matrix_w1_R6");
}
