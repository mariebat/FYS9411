#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <armadillo>
using namespace arma;


class MonteCarlo
{
public:
    MonteCarlo();
    MonteCarlo(int nEl, bool z, double A, int ncycles, double Omega, double ts); //use for importance sampling
    MonteCarlo(int nEl, bool z, double A, int ncycles, double Omega, double ts, double alph, double bet); //use for importance sampling


    void runMonteCarloInt(); //regular VMC integration program
    void MonteCarloSampling(int NumberMCsamples, double &cumulative_e, double &cumulative_e2); //VMC integration, optimized to run in parallel

    vec steepestDescent(vec x0, double gamma); //steepest descent


private:
    //Programs for steepest descent
    vec runMonteCarloInt_SD(vec guess); //returns expectation values of energy and wavefunciton
    vec gradLocalEnergy(vec expValues); // get derivative of local energy, w.r.t. alpha and beta
    double gradPsiAlpha(const mat &r); //get gradient of psi w.r.t. alpha
    double gradPsiBeta(const mat &r); //get gradient of psi w.r.t. beta
    double SPpsiGradAlpha(int index, double rx, double ry); //single partivcle gradient w.r.t. alpha

    //Local energy at specific position
    vec computeLocalEnergy_num(const mat &r);  //numeric derivative
    double computeLocalEnergy(const mat &r); //analytic expression, returns only local energy
    vec computeLocalEnergyvec(const mat &r); //analytic expression, returns vector with local energy, Ekin and Epot

    //wavefunction, at specific position
    double wavefunction(mat r); //trial wavefunction
    double SPwavefunction(double rx, double ry, int i); //single particle wavefunction
    double calcHermite(int n, double Z); //Hermite polynomial

    //helping functions
    double RelativeDistance(mat r, int i, int j);
    double singleparticle_posSQ(mat r, int i);
    int factorial(int A);

    //Quantum force, numerical derivative first and analytical expression second. Updates input QF matrix
    void quantumForce_num(const mat &r, mat &QF); //using numeric derivative
    void quantumForce(const mat &r, mat &QF); //using analytic expressions

    //Functions to find derivative
    double SPpsiDer(int index, double rx, double ry, int k); //derivative of Single particle wavefunction
    double SPpsiDoubleDer(int index, mat r); //double derivative of SP wavefunction
    double SlaterLaplace(int orbitals, mat r); //laplacian of slater determinant
    double hermiteDer(int n, double Z); //first derivative of hermite polynomial
    double hermiteDoubleDer(int n, double Z); //double derivative of hermite polynomial
    double JastrowDerivative(mat r, int i, int j, int k); //first derivative of jastrow factor
    double JastrowLaplacian(mat r); //laplacian of jastrow factor

    //System-specific variables
    int nDims;
    int nElectrons;
    double omega;

    double a; //Jastrow factor, turn off or on
    bool Z; //Coulomb interaction (to turn it off)

    //MC cycle parameters
    double stepLength; //for positions in brute force MC cycle
    double h; //For numerical derivative in local energy and quantum force
    double h2; // 1/(h*h)
    int nCycles;
    double timestep; //For updating position with importance sampling - parameter in soln to Langevin eqn

    //Variational parameters
    double alpha;
    double beta;

    //Position of electrons
    mat rOld;
    mat rNew;

    //Quantum force
    mat oldQF;
    mat newQF;

    double D; //Diffusion constant - equal to 1/2 in a.u.

    int nx; //spatial quantum number
    int ny; //spatial quantum number
    int sigma; //spin quantum number

    void map(int index); //maps orbital index to quantum numbers

    double pi = 3.14159265359;

};

#endif // MONTECARLO_H
