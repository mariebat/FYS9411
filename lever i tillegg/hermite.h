#ifndef HERMITE_H
#define HERMITE_H

//This class was initially created with the goal of computing the Hermite polynomials 
//for the Cartesian harmonic oscillator basis. As the scope of the project was altered, 
//the class was not completed. As of now, it only contains a function designed to compute 
//the Hermite polynomial at a given spatial coordinate Z (sqrt(w)*x or sqrt(w)*y) and 
//quantum number nx or ny. 

class hermite
{
public:
    hermite();
	hermite(int n, double Z);

    double calcHermite(int n,double Z);

private: 
	double H; 
};

#endif // HERMITE_H
