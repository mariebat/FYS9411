#include "hermite.h"
#include "harmonicoscillator.h"

#ifndef GAUSSHERMITE_H
#define GAUSSHERMITE_H


// This class consists of OpenMP code for numerical integration with Gauss-Hermite quadrature. 
//Here, the Coulomb integral matrix elements are computed. 

//The integrals are computed in a Cartesian basis. 

//The code is grouped into a class so it can be easily included in or removed from a Hartree Fock 
//or Monte Carlo calculated. 

//The gaussHermite class uses the Hermite polynomial class, as well as the harmonic oscillator class. 
//It needs access to both the qunatum number mapping functions of the harmonic oscillator class, as
//well as to the Hermite polynomials in the hermite class. 

class gaussHermite
{
public:
    gaussHermite();
    gaussHermite(double);

    // The integrand in Cartesian coordinates
    double Integrand(double, double,double, double, int, int, int, int);
    // The Gauss Hermite integration function
    double  GaussHermiteIntegration(int, int, int, int, int);

    // Getting the Gaussian quadrature weights and integration points
    void GaussHermiteQuadrature(double *, double *, int);

    void GaussLaguerreQuadrature(double *, double *, int, double);

    void GaussLegendreQuadrature(double, double, double *, double *, int);

    double GammaFunction( double);

    int factorial(int n);

    double int_function(double x);

private:
    hermite H;
    harmonicOscillator harm;
    double omega;
    double pi = 3.14159265359;

};


#endif // GAUSSHERMITE_H
