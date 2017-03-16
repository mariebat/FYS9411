#ifndef COULOMBINT_H
#define COULOMBINT_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

//#include <math.h>
//#include <string>
//#include <cstring>
//#include <sstream>
//#include <cstdarg>
//#include <algorithm>


//This class was made from the programs supplied by Morten Hjorth-Jensen. 
//The functions in the class generate the two-body interaction matrix elements. 

//The class uses an analytical expression for the Coulomb integral in 
//polar coordinates. 

//The Coulomb programs could have been directly implemented into the Hartree-Fock class, 
//but were separated in case there was a need to swiftly switch between Cartesian and polar corrdinates. 

class coulombInt
{
public:

    double Coulomb_HO(double &hw, int &ni, int &mi, int &nj, int &mj, int &nk, int &mk, int &nl, int &ml);
    double logratio1(int &int1, int &int2, int &int3, int &int4);
    double logratio2(int &G);
    double logfac(int &n);
    double product1(int &n1, int &m1, int &n2, int &m2, int &n3, int &m3, int &n4, int &m4);
    double logproduct2(int &n1, int &m1, int &n2, int &m2, int &n3, int &m3, int &n4, int &m4, int &j1, int &j2, int &j3, int &j4);
    double logproduct3(int &l1, int &l2, int &l3, int &l4, int &g1, int &g2, int &g3, int &g4);

private:


};


#endif // COULOMBINT_H
