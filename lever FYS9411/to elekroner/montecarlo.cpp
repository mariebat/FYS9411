#include "montecarlo.h"

#include <random>
#include <iomanip>
#include <math.h>

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

MonteCarlo::MonteCarlo()
{
    nDims = 2;
    alpha = 0.9903319732;
    //alpha = 0.68; // no jastrow optimal
    beta = 0.3951706571;

    stepLength = 1.5;
    nElectrons = 2;
    omega = 1.0;
    h = 0.001;
    h2 = 1/(h*h);
    timestep = 0.05; //0.001

    nCycles = 1e6;
    a = 1.0; //With Jastrow, two-electron case
    Z = 1; //With Coulomb interaction, two-electron QD

    D = 0.5;
}

MonteCarlo::MonteCarlo(int nEl, double steplength, int z, double A, int ncycles, double Omega)
{
    //Standard variables
    nDims = 2;
    alpha = 0.9903319732;
    //alpha = 0.68; // no jastrow optimal
    beta = 0.3951706571;

    //For numerical derivation in local energy
    h = 0.001;
    h2 = 1/(h*h);

    //Variables determined by the user
    stepLength = steplength;
    nElectrons = nEl;
    omega = Omega;
    timestep = 0.05; //0.001

    nCycles = ncycles;
    a = A; //with Jastrow?
    Z = z; //with Coulomb?

    D = 0.5;//diffusion constant
}

MonteCarlo::MonteCarlo(int nEl, int z, double A, int ncycles, double Omega, double ts)
{
    //Standard variables
    nDims = 2;
    alpha = 0.9903319732;
    //alpha = 0.68; // no jastrow optimal
    beta = 0.3951706571;

    stepLength = 1.5;

    //For numerical derivation in local energy
    h = 0.001;
    h2 = 1/(h*h);

    //Variables determined by the user
    nElectrons = nEl;
    omega = Omega;
    timestep = ts;

    nCycles = ncycles;
    a = A; //with Jastrow? 0 is no, 1 is yes
    Z = z; //with Coulomb? 0 is no, 1 is yes

    D = 0.5;//diffusion constant
}

MonteCarlo::MonteCarlo(int nEl, int z, double A, int ncycles, double Omega, double ts, double alph, double bet)
{
    //Standard variables
    nDims = 2;
    alpha = alph;
    beta = bet;

    stepLength = 1.5;

    //For numerical derivation in local energy
    h = 0.001;
    h2 = 1/(h*h);

    //Variables determined by the user
    nElectrons = nEl;
    omega = Omega;
    timestep = ts;

    nCycles = ncycles;
    a = A; //with Jastrow? 0 is no, 1 is yes
    Z = z; //with Coulomb? 0 is no, 1 is yes

    D = 0.5;//diffusion constant
}

void MonteCarlo::runMonteCarloInt_bruteForce()
{
    //Function to call Monte Carlo calculations.
    //Calls private functions comupteLocalEnergy and wavefunction, using the brute force Metropolis algorithm
    //no importance sampling!

    //Initialize random number generator
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);


    //Initialize position matrices to contain zeros
    rOld = zeros<mat>(nElectrons, nDims);
    rNew = zeros<mat>(nElectrons, nDims);

    //Initialize wavefunction and energies to zero
    double oldPsi = 0;
    double newPsi = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE;
    double acceptanceRate = 0;

    // initial trial positions
    for (int i = 0; i < nElectrons; i++) {
        for (int j = 0; j < nDims; j++) {
            rOld(i,j) = stepLength * (dist(mt) - 0.5);
        }
    }
    rNew = rOld;

    int runn = 0;
    //loop over MC cycles (actually run calc). Here: call private functions comupteLocalEnergy and wavefunction
    for (int run = 0; run < nCycles; run++)
    {
        oldPsi = wavefunction(rOld); //Store curent wavefunction value (double)

        //Test new position
        for (int i = 0; i < nElectrons; i++) {

            for (int j = 0; j < nDims; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(dist(mt) - 0.5);
            }
            newPsi = wavefunction(rNew);

            //Accept step? Yes--> update position, No --> reset position
            if (dist(mt) <= (newPsi*newPsi)/(oldPsi*oldPsi) ) {
                oldPsi = newPsi;
                for (int j = 0; j < nDims; j++) {
                    rOld(i,j) = rNew(i,j);
                }
                acceptanceRate++;

            }
            else{
                for (int j = 0; j < nDims; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }

            //Update the energies
            deltaE = computeLocalEnergy(rNew); //analytical expression for two electrons
            //deltaE = computeLocalEnergy_num(rNew); //numerical derivative to get local energy
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        } //end of loop over particles
        runn ++;


    }

    //Calculate energies and output results
    double energy = energySum/(nCycles * nElectrons);
    double energySquared = energySquaredSum/(nCycles * nElectrons);


    cout << setprecision(10) << "Energy: " << energy << ". " << endl;
    cout << setprecision(10) << "Energy (squared sum): " << energySquared << ". " << endl;

    //Compute and output variance
    double variance = energySquared - (energy*energy);
    cout << "Variance: " << variance << ". " << endl;

    //Output other relevant details
    cout << "Cycles: " << runn << ". " << endl;
    cout << "AcceptanceRate:" << acceptanceRate/(nCycles * nElectrons) << endl;
}

void MonteCarlo::runMonteCarloInt()
{
    //Function to call Monte Carlo calculations.
    //Calls private functions comupteLocalEnergy and wavefunction
    //using importance sampling

    //Initialize random number generator according to Guassian distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> gaussDist(0.0,1.0);
    std::uniform_real_distribution<double> uniDist(0.0, 1.0);

    //Store local energy
    // Open file for writing energies (for blocking) or positions (for one-body density)
    ofstream outfile;
    outfile.open("../pc/Dokumenter/MATLAB/FYS9411/project2/energies_E7_w1.txt", ios::out | ios::binary);
    outfile << setprecision(10);


    //Initialize position matrices to contain zeros
    rOld = zeros<mat>(nElectrons, nDims);
    rNew = zeros<mat>(nElectrons, nDims);

    //Initialize quantum force
    oldQF = zeros<mat>(nElectrons, nDims);
    newQF = zeros<mat>(nElectrons, nDims);

    vec Evec = zeros<vec>(4);

    //Initialize wavefunction and energies to zero
    double oldPsi = 0;
    double newPsi = 0;

    double energySum = 0;
    double energySquaredSum = 0;
    double r12 = 0;
    double Ekin = 0.0;
    double Epot = 0.0;

    double acceptanceRate = 0;

    // initial trial positions
    for (int i = 0; i < nElectrons; i++) {
        for (int j = 0; j < nDims; j++) {
            rOld(i,j) = gaussDist(gen)*sqrt(timestep);
        }
    }
    rNew = rOld;

    int runn = 0;
    //loop over MC cycles (actually run calc). Here: call private functions comupteLocalEnergy and wavefunction AND quantumForce
    for (int run = 0; run < nCycles; run++)
    {
        oldPsi = wavefunction(rOld); //Store curent wavefunction value (double)
        quantumForce(rOld,oldQF); //find quantum force, using analytical expressions
        //quantumForce_num(rOld,oldQF);

        //Set new position, one electron at a time, and test it
        for (int i = 0; i < nElectrons; i++) {
            for (int j = 0; j < nDims; j++) {
                rNew(i,j) = rOld(i,j) + gaussDist(gen)*sqrt(timestep) + oldQF(i,j)*timestep*D;
            }

            //Set position to old pos for the other particles; only move one particle at a time!
            for (int k = 0; k < nElectrons; k++) {
                if ( k != i) {
                    for (int j = 0; j < nDims; j++ ) {
                        rNew(k,j) = rOld(k,j);
                    }
                }
            }
            //Recalculate wavefunction and quantum force
            newPsi = wavefunction(rNew);
            quantumForce(rNew,newQF);
            //quantumForce_num(rNew,newQF);

            //Ratio of Greens functions for Metropolis-Hastings algorithm - numerical
            greenRatio = 0.0;
            for (int j = 0; j < nDims; j++) {
                greenRatio +=  0.5*(oldQF(i,j)+newQF(i,j))*( D*timestep*0.5*(oldQF(i,j)-newQF(i,j))-rNew(i,j)+rOld(i,j) );
            }

            //greenRatio = exp(0.5*greenRatio);
            greenRatio = exp(greenRatio);

            //Perform metropolis test, moving one particle at a time
            //Accept step? Yes--> update position and QF, No --> reset position and QF
            if (uniDist(gen) <= greenRatio*(newPsi*newPsi)/(oldPsi*oldPsi) ) {
                for (int j = 0; j < nDims; j++) {
                    rOld(i,j) = rNew(i,j);
                    oldQF(i,j) = newQF(i,j);
                }
                oldPsi = newPsi;
                acceptanceRate++;

            }
            else{
                for (int j = 0; j < nDims; j++) {
                    rNew(i,j) = rOld(i,j);
                    newQF(i,j) = oldQF(i,j);
                }
            }


            //Update the energies
            Evec =  computeLocalEnergyvec(rNew); //analytical expression for two electrons
            outfile << Evec(0) <<  "\t" << endl;

            energySum += Evec(0);
            energySquaredSum += Evec(0)*Evec(0);
            Ekin += Evec(1);
            Epot += Evec(2);
            r12 += Evec(3);

        } //end of loop over particles

        //outfile << r(0,0) << "\t" << r(0,1) << "\t" << r(1,0) << "\t" << r(1,1) << "\t" << endl;

        runn ++;
    }

    outfile.close();

    //Calculate energies and output results
    int n = nCycles*nElectrons;
    double energy = energySum/(n);
    double energySquared = energySquaredSum/(n);
    double meanDistance = r12/n;
    double K = Ekin/n;
    double V = Epot/n;

    cout << setprecision(10) << "Energy: " << energy << ". " << endl;
    cout << setprecision(10) << "Energy (squared sum): " << energySquared << ". " << endl;
    cout << setprecision(10) << "Mean distance between electrons: " << meanDistance << ". " << endl;
    cout << setprecision(10) << "Kinetic energy: " << K << ". " << endl;
    cout << setprecision(10) << "Potential energy: " << V << ". " << endl;

    //Compute and output variance
    double variance = energySquared - (energy*energy);
    double std = sqrt(variance/(n-1.0));
    cout << "Mean: " << energy << ". " << endl;
    cout << "Variance: " << variance << ". " << endl;
    cout << "Standard deviation: " << std << ". " << endl;

    //Output other relevant details
    cout << "Cycles: " << runn << ". " << endl;
    cout << "AcceptanceRate:" << acceptanceRate/(n) << endl;

}

vec MonteCarlo::steepestDescent(vec guess, double gamma)
{
    //function to perform steepest descent

    int maxDescents = 30;
    int i = 0;
    int dim = guess.n_elem;  //2!
    cout << "Starting steepest descent to optimize alpha, beta. " << endl;
    cout << "Dimensions: " << dim << endl;

    const double tolerance = 1.0e-6;

    vec gradient = zeros<vec>(2);
    vec expvec = zeros<vec>(5);

    vec x = zeros<vec>(dim);  //will contain alpha and beta

    double diffA;
    double diffB;

    while (i<=maxDescents ) //(i <= maxDescents)
    {


        expvec = runMonteCarloInt_SD(guess); //calculates local energy and gradients of wavefunction w.r.t. alpha and beta
        gradient = gradLocalEnergy(expvec); //calculates gradient of local energy w.r.t. alpha and beta

        x = guess - gamma*gradient;

        diffA = abs(x(0)-guess(0));
        diffB = abs(x(1)-guess(1));

        if ( diffA < tolerance && diffB < tolerance ) {
            break;
        }

        guess = x;

        i++;

        cout << "Just finished cycle #: " << i << endl;
    }

    cout << "DiffA: " << diffA << endl;
    cout << "DiffB: " << diffB << endl;

    return x; //contains (hopefully) optimal gamma and beta

}

vec MonteCarlo::gradLocalEnergy(vec expValues)
{
    //calculates gradient of local energy w.r.t alpha and beta
    //needs local energy and gradient of wavefunction w.r.t. alpha and beta

    vec ret = zeros<vec>(2);

    double EL = expValues(0);
    double dPsidA = expValues(1);
    double dPsidB = expValues(2);
    double prodA = expValues(3);
    double prodB = expValues(4);

    double gradA = 2*(prodA - dPsidA*EL);
    double gradB = 2*(prodB - dPsidB*EL);

    ret(0) = gradA;
    ret(1) = gradB;

    return ret;
}

vec MonteCarlo::runMonteCarloInt_SD(vec guess)
{
    //Program to run MC cycles within the steepest descent method
    //The program returns the expectation value of the energy

    //Initialize random number generator according to Guassian distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> gaussDist(0.0,1.0);
    std::uniform_real_distribution<double> uniDist(0.0, 1.0);


    //allocate room for vector to contain:
    //expectation value of local energy, expectation value of the derivative w.r.t alpha,
    //and exp. value of EL w.r.t. beta
    vec expvec = zeros<vec>(5);

    //use accepted vector to set private variables alpha, beta for each cycle
    alpha = guess(0);
    beta = guess(1);

    //Initialize position matrices to contain zeros
    rOld = zeros<mat>(nElectrons, nDims);
    rNew = zeros<mat>(nElectrons, nDims);

    //Initialize quantum force
    oldQF = zeros<mat>(nElectrons, nDims);
    newQF = zeros<mat>(nElectrons, nDims);

    //Initialize wavefunction and energies to zero
    double oldPsi = 0;
    double newPsi = 0;

    double energySum = 0;
    double derSumA = 0;
    double derSumB = 0;
    double prodSumA = 0;
    double prodSumB = 0;
    double deltaE, dPsidA, dPsidB;

    //double covar = 0.0;

    double acceptanceRate = 0;

    // initial trial positions
    for (int i = 0; i < nElectrons; i++) {
        for (int j = 0; j < nDims; j++) {
            rOld(i,j) = gaussDist(gen)*sqrt(timestep);
        }
    }
    rNew = rOld;

    int runn = 0;
    //loop over MC cycles (actually run calc). Here: call private functions comupteLocalEnergy and wavefunction AND quantumForce
    for (int run = 0; run < nCycles; run++)
    {
        oldPsi = wavefunction(rOld); //Store curent wavefunction value (double)
        quantumForce(rOld,oldQF); //find quantum force, using analytical expressions
        //quantumForce_num(rOld,oldQF);

        //Set new position, one electron at a time, and test it
        for (int i = 0; i < nElectrons; i++) {
            for (int j = 0; j < nDims; j++) {
                rNew(i,j) = rOld(i,j) + gaussDist(gen)*sqrt(timestep) + oldQF(i,j)*timestep*D;
            }

            //Set position to old pos for the other particles; only move one particle at a time!
            for (int k = 0; k < nElectrons; k++) {
                if ( k != i) {
                    for (int j = 0; j < nDims; j++ ) {
                        rNew(k,j) = rOld(k,j);
                    }
                }
            }
            //Recalculate wavefunction and quantum force
            newPsi = wavefunction(rNew);
            quantumForce(rNew,newQF);
            //quantumForce_num(rNew,newQF);

            //Ratio of Greens functions for Metropolis-Hastings algorithm - numerical
            greenRatio = 0.0;
            for (int j = 0; j < nDims; j++) {
                greenRatio +=  0.5*(oldQF(i,j)+newQF(i,j))*( D*timestep*0.5*(oldQF(i,j)-newQF(i,j))-rNew(i,j)+rOld(i,j) );
            }

            //greenRatio = exp(0.5*greenRatio);
            greenRatio = exp(greenRatio);

            //Perform metropolis test, moving one particle at a time
            //Accept step? Yes--> update position and QF, No --> reset position and QF
            if (uniDist(gen) <= greenRatio*(newPsi*newPsi)/(oldPsi*oldPsi) ) {
                for (int j = 0; j < nDims; j++) {
                    rOld(i,j) = rNew(i,j);
                    oldQF(i,j) = newQF(i,j);
                }
                oldPsi = newPsi;
                acceptanceRate++;

            }
            else{
                for (int j = 0; j < nDims; j++) {
                    rNew(i,j) = rOld(i,j);
                    newQF(i,j) = oldQF(i,j);
                }
            }

            //Update the energies
            deltaE = computeLocalEnergy(rNew); //analytical expression for two electrons
            dPsidA = gradPsiAlpha(rNew); //get gradient of wavefunction w.r.t. alpha (analytic)
            dPsidB = gradPsiBeta(rNew); //get gradient of wavefunction w.r.t. beta (analytic)

            energySum += deltaE;
            derSumA += dPsidA;
            derSumB += dPsidB;
            prodSumA += dPsidA*deltaE;
            prodSumB += dPsidB*deltaE;
        }

        runn ++;
    }


    //Calculate expectation value of energy for given alpha, beta and
    //return result
    int n = nCycles*nElectrons;

    double energy = energySum/(n); //expectation value of energy
    double dAlpha = derSumA/n;
    double dBeta = derSumB/n;
    double prodA = prodSumA/n;
    double prodB = prodSumB/n;

    expvec(0) = energy;
    expvec(1) = dAlpha;
    expvec(2) = dBeta;
    expvec(3) = prodA;
    expvec(4) = prodB;

    return expvec;
}

double MonteCarlo::gradPsiAlpha(const mat &r)
{
    //calculates gradient of wavefunction w.r.t. alpha from analytic expression
    double sum = 0.0;


    for (int i = 0; i < nElectrons; i++) {
        double A = 0.0;
        for (int j = 0; j < nDims; j++) {
            A += r(i,j) * r(i,j);
        }
        sum += A;
    }

    double gradAlpha = -0.5*omega*sum;
    return gradAlpha;
}

double MonteCarlo::gradPsiBeta(const mat &r)
{
    //calculates gradient of wavefunction w.r.t. beta form analytic expression
    double r12 = 0.0;

    for (int i = 0; i < nElectrons; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            double u = 0.0;
            for (int k = 0; k < nDims; k++) {
                u += (r(i,k) - r(j,k) ) * (r(i,k) - r(j,k) );
            }
            r12 += sqrt(u);
        }
    }

    double gradBeta = -a*r12*r12/((1+beta*r12)*(1+beta*r12));
    return gradBeta;
}

void MonteCarlo::quantumForce(const mat &r, mat &QF) //analytic  expression
{
    //calculates quantum force from analytic expression as fn. of position
    double r1r2 = 0.0;
    double u;

    for (int i = 0; i < nElectrons; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            u = 0.0;
            for (int k = 0; k < nDims; k++) {
                u += (r(i,k) - r(j,k) ) * (r(i,k) - r(j,k) );
            }
            r1r2 +=  sqrt(u);
        }
    }

    for (int i = 0; i < nElectrons; i++) {
        for (int j = nElectrons-1; j > -1; j--) {
            for (int k = 0; k < nDims; k++) {
                QF(i,k) = -2*alpha*omega*r(i,k) + 2*a*((r(i,k) - r(j,k))/(r1r2*(1+beta*r1r2)*(1+beta*r1r2) ));
            }
        }
    }

}

void MonteCarlo::quantumForce_num(const mat &r, mat &QF) //numerical derivative
{
    //calculates quantum force using numerical derivative

    //Initialize
    mat rPlus = zeros<mat>(nElectrons, nDims);
    mat rMin = zeros<mat>(nElectrons, nDims);

    rPlus = rMin = r;

    double PsiMin = 0.0;
    double PsiPlus = 0.0;

    double PsiNow = wavefunction(r);

    for(int i = 0; i < nElectrons; i++) {
        for (int j = 0; j < nDims; j++) {
            rPlus(i,j) += h;
            rMin(i,j) -= h;
            PsiPlus = wavefunction(rPlus);
            PsiMin = wavefunction(rMin);

            QF(i,j) = (PsiPlus - PsiMin);

            rPlus(i,j) = r(i,j);
            rMin(i,j) = r(i,j);

        }
    }

    QF = 2*h*QF/(2*PsiNow);


}

double MonteCarlo::wavefunction(mat r)
{
    //calculates the wavefunction at specific position r

    double C = 1.0; //normalization constant
    double first = 0.0;
    double u = 0.0;
    double r1r2 = 0.0;


    //First exponential
    for (int i = 0; i < nElectrons; i++) {
        double rSingleEl = 0.0;
        for (int j = 0; j < nDims; j++) {
            rSingleEl += r(i,j) * r(i,j);
        }
        first += rSingleEl;
    }

    //Second exponential
    for (int i = 0; i < nElectrons; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            r1r2 = 0.0;
            for (int k = 0; k < nDims; k++) {
                r1r2 += (r(i,k) - r(j,k) ) * (r(i,k) - r(j,k) );
            }
            u +=  sqrt(r1r2);
        }
    }
    //double rij = sqrt(u);
    double rij = u;

    double Jastrow = a*rij/(1 + (beta*rij));

    //Wavefunction; Can turn of Jastrow factor by setting a = 0 at beginning. Unperturbed: alpha = 1
    return C*exp(-0.5*alpha*omega*first)*exp(Jastrow);
}

vec MonteCarlo::computeLocalEnergyvec(const mat &r)
{
    //calculates local energy from analytic expression
    //Stores local energy, kinetic energy, potential energy and mean electron-electron distance in vector

    //Initialize
    double EKin;
    double EPot;
    double Eelel = 0.0;

    double u = 0.0;
    double r1r2 = 0.0;

    //Potential energy - harmonic oscillator potential! mulitply with 1/2 omega^2 at the end
    EPot =  0.5*omega*omega*accu(r%r);

    //Contribution from electron-electron potential (repulsion) - also to kinetic energy
    for (int i = 0; i < nElectrons; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            r1r2 = 0.0;
            for (int k = 0; k < nDims; k++) {
                r1r2 += (r(i,k) - r(j,k) ) * (r(i,k) - r(j,k) );
            }
            u +=  sqrt(r1r2);
        }
    }
    //double rij = sqrt(u);
    double rij = u;
    Eelel = Z / rij; //Can set Z = 0 to turn off Coulomb repulsion. Otherwise, it is 1


    //Kinetic energy (analytical expression, two electrons)
    double nor12 = alpha*alpha*omega*omega*accu(r%r) - 4*alpha*omega;
    double withr12 = (2*a/((1+beta*rij)*(1+beta*rij)))*( (a/((1+beta*rij)*(1+beta*rij))) + (1/rij) - (2*beta/(1+beta*rij)) );

    EKin = nor12 - (2*a*alpha*omega*rij/( (1+beta*rij)*(1+beta*rij)) ) + withr12;

    double EL = (Eelel - 0.5*EKin + EPot );
    vec Evec = zeros<vec>(4);

    Evec(0) = EL;
    Evec(1) = -0.5*EKin;
    Evec(2) = EPot+Eelel;
    Evec(3) = rij;

    return Evec;
}

double MonteCarlo::computeLocalEnergy(const mat &r)
{
    //Compute local enegy from analytical expression

    //Initialize
    double EKin;
    double EPot;
    double Eelel = 0.0;

    double U = accu(r%r);
    double u = 0.0;
    double r1r2 = 0.0;

    //Potential energy - harmonic oscillator potential! mulitply with 1/2 omega^2 at the end
    EPot =  0.5*omega*omega*U;

    //Contribution from electron-electron potential (repulsion) - also to kinetic energy
    for (int i = 0; i < nElectrons; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            r1r2 = 0.0;
            for (int k = 0; k < nDims; k++) {
                r1r2 += (r(i,k) - r(j,k) ) * (r(i,k) - r(j,k) );
            }
            u +=  sqrt(r1r2);
        }
    }
    //double rij = sqrt(u);
    double rij = u;
    Eelel = Z / rij; //Can set Z = 0 to turn off Coulomb repulsion. Otherwise, it is 1


    //Kinetic energy (analytical expression, two electrons)
    double nor12 = alpha*alpha*omega*omega*U - 4*alpha*omega;
    double withr12 = (2*a/((1+beta*rij)*(1+beta*rij)))*( (a/((1+beta*rij)*(1+beta*rij))) + (1/rij) - (2*beta/(1+beta*rij)) );

    EKin = nor12 - (2*a*alpha*omega*rij/( (1+beta*rij)*(1+beta*rij)) ) + withr12;

    return (Eelel - 0.5*EKin + EPot );
}

double MonteCarlo::computeLocalEnergy_num(const mat &r)
{
    //Compute local energy using numerical derivatives

    //Initialize
    mat rPlus = zeros<mat>(nElectrons, nDims);
    mat rMin = zeros<mat>(nElectrons, nDims);
    rPlus = r;
    rMin = r;

    double psiMin = 0.0;
    double psiPlus = 0.0;
    double psiNow = wavefunction(r);

    double EKin = 0.0;
    double EPot;
    double Eelel = 0.0;

    double u = 0.0;
    double r1r2 = 0.0;
    double r1 = 0.0;
    double U = 0.0;

    //Potential energy - harmonic oscillator potential! mulitply with 1/2 omega^2 at the end
    for (int i = 0; i < nElectrons; i++) {
        r1 = 0.0;
        for (int j = 0; j < nDims; j++) {
            r1 += r(i,j)*r(i,j);
        }
        U += r1;
    }
    EPot =  0.5*omega*omega*U;

    //Contribution from electron-electron potential (repulsion)
    for (int i = 0; i < nElectrons; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            r1r2 = 0.0;
            for (int k = 0; k < nDims; k++) {
                r1r2 += (r(i,k) - r(j,k) ) * (r(i,k) - r(j,k) );
            }
            u += sqrt(r1r2);
        }
    }
    double rij = u;
    Eelel = Z / rij; //Can set Z = 0 to turn off Coulomb repulsion. Otherwise, it is 1

    //Kinetic energy (numerical derivative - brute force! )
    for (int i = 0; i < nElectrons; i++) {
        for (int j = 0; j < nDims; j++) {
            rPlus(i,j) += h;
            rMin(i,j) -= h;

            psiMin = wavefunction(rMin);
            psiPlus = wavefunction(rPlus);

            EKin -= (psiMin + psiPlus - 2*psiNow);

            rPlus(i,j) = r(i,j);
            rMin(i,j) = r(i,j);
        }
    }

    EKin = 0.5 * h2 * EKin / psiNow;

    return (EPot + EKin + Eelel);
}




