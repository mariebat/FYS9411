#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <armadillo>
using namespace arma;


class MonteCarlo
{
public:
    MonteCarlo();
    MonteCarlo(int nEl, double steplength, int z, double A, int ncycles, double Omega); //use for brute force MC
    MonteCarlo(int nEl, int z, double A, int ncycles, double Omega, double ts); //use for importance sampling
    MonteCarlo(int nEl, int z, double A, int ncycles, double Omega, double ts, double alph, double bet); //use for importance sampling

    void runMonteCarloInt_bruteForce(); //brute force Metropolis algorithm
    void runMonteCarloInt(); //Metropolis-Hastings algorithm with importance sampling

    vec steepestDescent(vec x0, double gamma); //program to run Steepest descent


private:
    //program sfor steepest descent
    vec runMonteCarloInt_SD(vec guess); //run MC within SD; returns expectation values of energy and wavefunciton
    vec gradLocalEnergy(vec expValues); // get derivative of local energy, w.r.t. alpha and beta

    double gradPsiAlpha(const mat &r); //get gradient of psi w.r.t. alpha
    double gradPsiBeta(const mat &r); //get gradient of psi w.r.t. beta

    //local energy, numeric derivative and analytic expression, respectively, at specific position
    double computeLocalEnergy_num(const mat &r);
    double computeLocalEnergy(const mat &r);
    vec computeLocalEnergyvec(const mat &r); //returns vector, with kinetic and potential energies

    //wavefunction, at specific position
    double wavefunction(mat r);

    //Quantum force, numerical derivative first and analytical expression second. Updates input QF matrix
    void quantumForce_num(const mat &r, mat &QF);
    void quantumForce(const mat &r, mat &QF);

    //System-specific variables
    int nDims;
    int nElectrons;
    double omega;

    double a; //Jastrow
    int Z; //Coulomb interaction (to turn it off)

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

    //Importance sampling parameters
    double D; //Diffusion constant - equal to 1/2 in a.u.
    double greenRatio; //Ratio of green's functions for importance sampling

};

#endif // MONTECARLO_H
