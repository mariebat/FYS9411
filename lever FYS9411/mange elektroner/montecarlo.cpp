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
    a = 1; //With Jastrow, two-electron case
    Z = 1; //With Coulomb interaction, two-electron QD

    D = 0.5;
}

MonteCarlo::MonteCarlo(int nEl, bool z, double A, int ncycles, double Omega, double ts)
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

MonteCarlo::MonteCarlo(int nEl,bool z, double A, int ncycles, double Omega, double ts, double alph, double bet)
{
    //Standard variables
    nDims = 2;
    D = 0.5; //diffusion constant
    stepLength = 1.5;

    //For numerical derivation in local energy
    h = 0.001;
    h2 = 1/(h*h);

    //Variables determined by the user
    nElectrons = nEl;
    omega = Omega;
    timestep = ts;

    alpha = alph;
    beta = bet;

    nCycles = ncycles;
    a = A; //with Jastrow? 0 is no, 1 is yes
    Z = z; //with Coulomb? 0 is no, 1 is yes

}

int MonteCarlo::factorial(int A)
{
  return (A == 1 || A == 0) ? 1 : factorial(A - 1) * A;
}

void MonteCarlo::MonteCarloSampling(int NumberMCsamples, double &cumulative_e, double &cumulative_e2)
{
    //Function to call Monte Carlo calculations.
    //Calls private functions comupteLocalEnergy and wavefunction
    //using importance sampling
    //For paralellized calculations

    //Initialize random number generator according to Guassian distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> gaussDist(0.0,1.0);
    std::uniform_real_distribution<double> uniDist(0.0, 1.0);

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
    double energySquaredSum = 0;
    double EL;


    // initial trial positions
    for (int i = 0; i < nElectrons; i++) {
        for (int j = 0; j < nDims; j++) {
            rOld(i,j) = gaussDist(gen)*sqrt(timestep);
        }
    }
    rNew = rOld;

    //loop over MC cycles (actually run calc). Here: call private functions comupteLocalEnergy and wavefunction AND quantumForce
    for (int run = 0; run < NumberMCsamples; run++)
    {

        if ( run%100000 == 0)
            cout << "Run number: " << run << endl;

        oldPsi = wavefunction(rOld); //Store curent wav efunction
        quantumForce(rOld,oldQF); //find quantum force, using analytical expressions
        //quantumForce_num(rOld,oldQF);
        //Set new position, one electron at a time, and    cout << "hit da" << endl; test it
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
            double greenRatio = 0.0;
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

            }

        } //end of loop over particles

        //Calculate relevant parameters;
        //local enegy, energy squared
        EL =  computeLocalEnergy(rOld); //analytical expression for local energy for two electrons

        energySum += EL;
        energySquaredSum += EL*EL;

    }
    cumulative_e = energySum/NumberMCsamples; //mean
    cumulative_e2 = energySquaredSum/NumberMCsamples; //mean
}

void MonteCarlo::runMonteCarloInt()
{
    //Function to call Monte Carlo calculations.
    //Calls private functions comupteLocalEnergy and wavefunction
    //using importance sampling

    //Initialize random number generator according to Guassian distribution
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::normal_distribution<double> gaussDist(0.0,1.0);
    std::uniform_real_distribution<double> uniDist(0.0,1.0);

    //Store local energy
    // Open file for writing, writing results in formated output for plotting.
    //Can write either energy, for blocking, or electron positions, for one-body density
    ofstream outfile;
    outfile.open("../../../../MATLAB/FYS9411/project2/energies_6el_w1.txt", ios::out | ios::binary);
    outfile << setprecision(10);


    //Initialize position matrices to contain zeros
    rOld = zeros<mat>(nElectrons, nDims);
    rNew = zeros<mat>(nElectrons, nDims);

    //Initialize quantum force
    oldQF = zeros<mat>(nElectrons, nDims);
    newQF = zeros<mat>(nElectrons, nDims);

    //Initialize energy sums to zero
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
    //rNew = rOld;
    double oldPsi = wavefunction(rOld);
    quantumForce(rOld,oldQF);

    //int runn = 0;
    //loop over MC cycles (actually run calc).
    //Here: call private functions comupteLocalEnergy and wavefunction AND quantumForce
    for (int run = 0; run < nCycles; run++)
    {

        if ( run%10000 == 0)
           cout << "Run number: " << run << endl;

        //oldPsi = wavefunction(rOld); //Store curent wav efunction
        //quantumForce(rOld,oldQF); //find quantum force, using analytical expressions

        //Set new position, one electron at a time
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
            double newPsi = wavefunction(rNew);
            quantumForce(rNew,newQF);

            //Ratio of Greens functions for Metropolis-Hastings algorithm - numerical
            double greenRatio = 0.0;
            for (int j = 0; j < nDims; j++) {
                greenRatio +=  0.5*(oldQF(i,j)+newQF(i,j))*( D*timestep*0.5*(oldQF(i,j)-newQF(i,j))-rNew(i,j)+rOld(i,j) );
            }

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

        } //end of loop over particles

        //Calculate relevant parameters;
        //local enegy, energy squared, kinetic energy, potential energy, average distance between electrons

        vec Evec =  computeLocalEnergyvec(rOld); //analytical expression for two electrons

        energySum += Evec(0);
        energySquaredSum += Evec(0)*Evec(0);
        Ekin += Evec(1);
        Epot += Evec(2);
        r12 += Evec(3);


        //Write EL to file to use for blocking
        outfile << Evec(0) << "\t" << endl;


        //Skriv posisjoner til fil for one-body density
//        for (int i = 0; i < nElectrons; i++) {
//            outfile << rOld(i,0) << ", " << rOld(i,1) << ", ";
//        }
//        outfile << endl;



    }

    outfile.close();

    //Calculate energies and output results
    int n = nCycles;//*nElectrons;
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
    cout << "Cycles: " << nCycles << ". " << endl;
    cout << "AcceptanceRate:" << acceptanceRate/(nElectrons*nCycles) << endl;




}

vec MonteCarlo::steepestDescent(vec guess, double gamma)
{
    //Program to run steepest descent calculations

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

        expvec = runMonteCarloInt_SD(guess); //run MC integration; get local energy and wavefunction grdients
        gradient = gradLocalEnergy(expvec); //gradient of local energy w.r.t. alpha and beta

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
    //Calculates gradient of local energy w.r.t. alpha and beta for steepest descent
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
    //Function which runs MC cycles for the SD method
    //The program returns the expectation value of the energy and gradients of the wavefunction
    //w.r.t. alpha and beta

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

    // initial trial positions
    for (int i = 0; i < nElectrons; i++) {
        for (int j = 0; j < nDims; j++) {
            rOld(i,j) = gaussDist(gen)*sqrt(timestep);
        }
    }
    oldPsi = wavefunction(rOld); //Store curent wavefunction value (double)
    quantumForce(rOld,oldQF); //find quantum force, using analytical expression

    //loop over MC cycles (actually run calc). Here: call private functions comupteLocalEnergy and wavefunction AND quantumForce
    for (int run = 0; run < nCycles; run++)
    {


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
            double greenRatio = 0.0;
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

            }


        }

        //Update the energies
        deltaE = computeLocalEnergy(rNew); //analytical expression for two electrons
        dPsidA = gradPsiAlpha(rNew);
        dPsidB = gradPsiBeta(rNew);

        energySum += deltaE;
        derSumA += dPsidA;
        derSumB += dPsidB;
        prodSumA += dPsidA*deltaE;
        prodSumB += dPsidB*deltaE;

    }


    //Calculate expectation value of energy for given alpha, beta and
    //return result
    int n = nCycles;

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
    //calculates gradient of wavefunction w.r.t alpha  from analytic expression

    int nOrbitals = nElectrons/2;

    mat SDup(nOrbitals,nOrbitals,fill::zeros);
    mat SDdown(nOrbitals,nOrbitals,fill::zeros);

    for (int i = 0; i < nOrbitals; i++) {
        for (int j = 0; j < nOrbitals; j++) {
            SDup(i,j) = SPwavefunction(r(i,0),r(i,1), 2*j);
            SDdown(i,j) = SPwavefunction(r(i+nOrbitals,0),r(i+nOrbitals,1), 2*j+1);
        }
    }

    mat SDupinv = inv(SDup);
    mat SDdowninv = inv(SDdown);


    // compute the gradient of the SD w.r.t. alpha
    double Ssum = 0.0;

    for (int i = 0; i < nElectrons; i++) {



        if (i < nOrbitals)
               for (int j = 0; j < nOrbitals; j++)
               {
                   Ssum += SPpsiGradAlpha(2*j,r(i,0),r(i,1))*SDupinv(j,i);
               }

        else
            for (int j = 0; j < nOrbitals; j++)
                {
                    Ssum += SPpsiGradAlpha(2*j+1,r(i,0),r(i,1))*SDdowninv(j,i-nOrbitals);
                }

    }

    return Ssum;


}

double MonteCarlo::SPpsiGradAlpha(int index, double rx, double ry)
{
    //calculates gradient of single-particle wavefunction w.r.t alpha

    map(index); //nx, ny

    double Hx, Hy, Hxd, Hyd;

    double u = sqrt(omega*alpha)*rx;
    double z = sqrt(omega*alpha)*ry;
    double k = -0.5*omega*alpha*(rx*rx + ry*ry);

    double du = sqrt(omega)*rx/(2*sqrt(alpha));
    double dz = sqrt(omega)*ry/(2*sqrt(alpha));
    double dk = -0.5*omega*(rx*rx + ry*ry);


    //Compute polynomials
    Hx = calcHermite(nx, u);
    Hy = calcHermite(ny, z);
    Hxd = hermiteDer(nx, u);
    Hyd = hermiteDer(ny, z);

    double derHxHy = Hxd*du*Hy + Hx*Hyd*dz;

    double dphi = derHxHy*exp(k) + Hx*Hy*dk*exp(k);

    return dphi;
}

double MonteCarlo::gradPsiBeta(const mat &r)
{
    //calculates gradient of wavefunction w.r.t beta from analytic expression

    double gradBeta = 0.0;
    double sum = 0.0;

    for (int i = 0; i < nElectrons-1; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            double rij = RelativeDistance(r,i,j);
            gradBeta -= a*rij*rij/((1.0+beta*rij)*(1.0+beta*rij));
            sum += a*rij/(1.0+beta*rij);
        }
    }

    return gradBeta;
}

double MonteCarlo::JastrowDerivative(mat r, int i, int j, int k)
{
    //calculates derivative of jastrow faxctor

    double rij = RelativeDistance(r, i, j);
    double derJ = (r(i,k) - r(j,k))*a/(rij*pow((1+beta*rij),2));

    return derJ;
}

double MonteCarlo::SPpsiDoubleDer(int index, mat r)
{
    //returns single particle double derivative in both dimensions

    map(index); //nx, ny

    double phider;
    double Hx, Hy, Hxd, Hyd, Hxdd, Hydd;

    double rx = r(index,0);
    double ry = r(index,1);

    //Copmpute amplitudes
    double Ax = 1; //pow((omega/pi),0.25)*(1/sqrt(pow(2,nx)*factorial(nx)));
    double Ay = 1; // pow((omega/pi),0.25)*(1/sqrt(pow(2,ny)*factorial(ny)));

    double pos2 = singleparticle_posSQ(r, index);

    //Compute polynomials
    Hx = calcHermite(nx, sqrt(omega*alpha)*rx);
    Hy = calcHermite(ny, sqrt(omega*alpha)*ry);
    Hxd = hermiteDer(nx, sqrt(omega*alpha)*rx);
    Hyd = hermiteDer(ny, sqrt(omega*alpha)*ry);
    Hxdd = hermiteDoubleDer(nx,sqrt(omega*alpha)*rx);
    Hydd = hermiteDoubleDer(ny,sqrt(omega*alpha)*ry);

    //Compute wavefunction
    double phiddx = Hy*(Hxdd+Hx*(omega*omega*alpha*alpha*rx*rx-omega*alpha)-2*omega*alpha*rx*Hxd);
    double phiddy = Hx*(Hydd+Hy*(omega*omega*alpha*alpha*ry*ry-omega*alpha)-2*omega*alpha*ry*Hyd);
    phider = Ax*Ay*(phiddx + phiddy)*exp(-0.5*alpha*omega*pos2);

    return phider;
}

double MonteCarlo::SPpsiDer(int index, double rx, double ry, int k)
{
    //returns derivative of Single particle waveunction in dimension k

    map(index); //nx, ny

    double phider = 0;
    double Hx, Hy, Hxd, Hyd;

    //Copmpute amplitudes
    double Ax = 1; //pow((omega/pi),0.25)*(1/sqrt(pow(2,nx)*factorial(nx)));
    double Ay = 1; // pow((omega/pi),0.25)*(1/sqrt(pow(2,ny)*factorial(ny)));

    //Compute polynomials
    Hx = calcHermite(nx, sqrt(omega*alpha)*rx);
    Hy = calcHermite(ny, sqrt(omega*alpha)*ry);

    double pos2 = rx*rx + ry*ry;
    if ( k == 0)
    {
        Hxd = hermiteDer(nx, sqrt(omega*alpha)*rx);
        phider = Ax*Ay*Hy*(Hxd*sqrt(alpha*omega) - Hx*alpha*omega*rx)*exp(-alpha*omega*0.5*pos2);

    }

    else if (k == 1)
    {
        Hyd = hermiteDer(ny, sqrt(omega*alpha)*ry);
        phider = Ax*Ay*Hx*(Hyd*sqrt(alpha*omega) - Hy*alpha*omega*ry)*exp(-alpha*omega*0.5*pos2);
    }


    else
    {
        cout << "Wrong k. Try again. " << endl;
    }

    return phider;

}

void MonteCarlo::quantumForce(const mat &r, mat &QF) //analytic  expression
{
    //calculates quantum force, updates QF

    int nOrbitals = nElectrons/2;

    mat SDup(nOrbitals,nOrbitals,fill::zeros);
    mat SDdown(nOrbitals,nOrbitals,fill::zeros);

    //fill A with Hermite polynomial wavefunctions
    for (int i = 0; i < nOrbitals; i++) {
        for (int j = 0; j < nOrbitals; j++) {
            SDup(i,j) = SPwavefunction(r(i,0),r(i,1), 2*j);
            SDdown(i,j) = SPwavefunction(r(i+nOrbitals,0),r(i+nOrbitals,1), 2*j+1);
        }
    }

    mat SDupinv = inv(SDup);
    mat SDdowninv = inv(SDdown);


    // compute the jastrow gradient
    for (int i = 0; i < nElectrons; i++) {
        for (int k = 0; k < nDims; k++) {

            // Jastrow factor contribution
            double Jsum = 0.0;
            double Ssum = 0.0;
            for (int j = 0; j < i; j++) {
                double rji = RelativeDistance(r,j,i);
                Jsum += ((r(i,k)-r(j,k))/(rji))*a/((1.0+beta*rji)*(1.0+beta*rji));
            }

            for (int j = i+1; j < nElectrons; j++) {
                double rij = RelativeDistance(r,i,j);
                Jsum -= ((r(j,k)-r(i,k))/(rij))*a/((1+beta*rij)*(1+beta*rij));
            }

//            for (int j = i+1; j < nElectrons; j++) {
//                double rij = RelativeDistance(r,i,j);
//                Jsum += (r(i,k) - r(j,k))*a/(rij*pow((1.0+beta*rij),2));
//            }

            if (i < nOrbitals)
                for (int j = 0; j < nOrbitals; j++)
                {
                    Ssum += SPpsiDer(2*j,r(i,0),r(i,1),k)*SDupinv(j,i);
                }

            else
                for (int j = 0; j < nOrbitals; j++)
                {
                    Ssum += SPpsiDer(2*j+1,r(i,0),r(i,1),k)*SDdowninv(j,i-nOrbitals);
                }
            QF(i,k) = 2.0*(Jsum + Ssum);

        }
    }

}

void MonteCarlo::quantumForce_num(const mat &r, mat &QF) //numerical derivative
{
    //calcultes quantum force from numerical derivative

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

double MonteCarlo::RelativeDistance(mat r, int i, int j)
{
    //calculates relative distancve between two electrons
    double r_ij = 0;
    for (int k = 0; k < nDims; k++) {
        r_ij += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
    }
    return sqrt(r_ij);
}

double MonteCarlo::singleparticle_posSQ(mat r, int i)
{
    //calcuates position squared of particle i

    double r_single_particle = 0;
    for (int j = 0; j < nDims; j++) {
      r_single_particle  += r(i,j)*r(i,j);
    }
    return r_single_particle;
}

double MonteCarlo::SPwavefunction(double rx, double ry, int i)
{
    //calculates single-particle wavefunction depending on position and orbital index

    map(i); //finds nx, ny and saves them


    //This function computes the Harmonic oscillator wavefunction phi at spatial coordinates x and y,
    //for specific quantum numbers nx and ny. The spin quantum number is not relevant. The harmonic oscillator frequency must be specified.

    double phi;
    double Hx;
    double Hy;

    //Copmpute amplitudes
    double Ax = 1; //  = pow((omega/pi),0.25)*(1/sqrt(pow(2,nx)*factorial(nx)));
    double Ay = 1; // pow((omega/pi),0.25)*(1/sqrt(pow(2,ny)*factorial(ny)));

    //Compute polynomials
    Hx = calcHermite(nx, sqrt(omega*alpha)*rx);
    Hy = calcHermite(ny, sqrt(omega*alpha)*ry);

    double pos2 = rx*rx + ry*ry;

    //Compute wavefunction
    phi = Ax*Ay*Hx*Hy*exp(-alpha*omega*0.5*pos2);

    return phi;
}

void MonteCarlo::map(int index)
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

double MonteCarlo::hermiteDer(int n, double Z)
{
    //Function to calculate the derivative of Hermite polynomials

    double H;
    switch (n)
    {
    case 0:
        H = 0;
        break;
    case 1:
        H = 2;
        break;
    case 2:
        H = 8 * Z;
        break;
    case 3:
        H = 8*3*Z*Z - 12;
        break;
    default:
        cout << "Wrong spatial quantum number given. It must be an integer between (or including) 0 and 2. Abort!";
        cout << "Returning -1 for Hermite polynomial derivatvive." << endl;
        cout << endl;
        return -1;
    }

    return H;

}

double MonteCarlo::hermiteDoubleDer(int n, double Z)
{
    //Function to calculate the double derivative of Hermite polynomials

    double H;
    switch (n)
    {
    case 0:
        H = 0;
        break;
    case 1:
        H = 0;
        break;
    case 2:
        H = 8;
        break;
    case 3:
        H = 8*3*2*Z;
        break;
    default:
        cout << "Wrong spatial quantum number given. It must be an integer between (or including) 0 and 2. Abort!";
        cout << "Returning -1 for Hermite polynomial double derivatvive." << endl;
        cout << endl;
        return -1;
    }

    return H;

}

double MonteCarlo::calcHermite(int n, double Z)
{
    //Function to calculate the Hermitian polynomial corresponding to the quantum number nx or ny,
    //both of which have values between (including) 0 and 9.

    // Z = sqrt(w)*x or sqrt(w)*y - anyway, a number, spatial coordinate
    double H;
    switch (n)
    {
    case 0:
        H = 1;
        break;
    case 1:
        H = 2 * Z;
        break;
    case 2:
        H = 4 * Z*Z - 2;
        break;
    case 3:
        H = 8 * Z*Z*Z - 12 * Z;
        break;
    case 4:
        H = 16 * pow(Z,4) - 48 * Z*Z + 12;
        break;
    case 5:
        H = 32 * pow(Z,5) - 160 * Z*Z + 120 * Z;
        break;
    case 6:
        H = 64 * pow(Z,6) - 480 * pow(Z,4) + 720 * Z*Z - 120;
        break;
    case 7:
        H = 128 * pow(Z,7) - 1344 * pow(Z,5) + 3360 * Z*Z*Z - 1680 * Z;
        break;
    case 8:
        H = 256 * pow(Z,8) - 3584 * pow(Z,6) + 13440 * pow(Z,4) - 13440 * Z*Z + 1680;
        break;
    case 9:
        H = 512 * pow(Z,9) - 9216 * pow(Z,7) + 48384 * pow(Z,5) - 80640 * Z*Z*Z + 30240 * Z;
        break;
    default:
        cout << "Wrong spatial quantum number given. It must be an integer between (or including) 0 and 2. Abort!";
        cout << "Returning -1 for Hermite polynomial." << endl;
        cout << endl;
        return -1;
    }

    return H;

}

double MonteCarlo::wavefunction(mat r)
{
    //returns trial wavefunction, including slater determinant and jastrow factor

    double wf = 0.0;
    int nOrbitals = nElectrons/2;
    mat SDup(nOrbitals,nOrbitals,fill::zeros);
    mat SDdown(nOrbitals,nOrbitals,fill::zeros);

    //fill A with Hermite polynomial wavefunctions
    //fill A with Hermite polynomial wavefunctions
    for (int i = 0; i < nOrbitals; i++) {
        for (int j = 0; j < nOrbitals; j++) {
            SDup(i,j) = SPwavefunction(r(i,0),r(i,1), 2*j);
            SDdown(i,j) = SPwavefunction(r(i+nOrbitals,0),r(i+nOrbitals,1), 2*j+1);
        }
    }
    wf = det(SDup)*det(SDdown);


    //Second exponential
    for (int i = 0; i < nElectrons-1; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            double rij = RelativeDistance(r, i, j);
            wf *= exp(a*rij/(1.0+beta*rij));
        }
    }

    return wf;
}

double MonteCarlo::SlaterLaplace(int orbitals, mat r)
{
    //returns laplacian of slater determinant

    double EHO;
    switch(orbitals)
    {
    case 1:
        EHO = 1;
        break;
    case 3:
        EHO = 5;
        break;
    case 6:
        EHO = 14;
        break;
    default:
        cout << "Wrong nOrbitals in SlaterLaplace - need closed shell systems." << endl;
        return -1;
    }
    double pos2 = accu(r%r);

    double slaterLaplace = alpha*alpha*omega*omega*pos2 - 4*alpha*omega*EHO;

    return slaterLaplace;
}

vec MonteCarlo::computeLocalEnergyvec(const mat &r)
{
    //returns local energy in vector with kinetic and potential energies
    //uses analytical expressions

    //Initialize
    int nOrbitals = nElectrons/2;

    mat SDup(nOrbitals,nOrbitals,fill::zeros);
    mat SDdown(nOrbitals,nOrbitals,fill::zeros);

    //Make Slater matrices
    for (int i = 0; i < nOrbitals; i++) {
        for (int j = 0; j < nOrbitals; j++) {
            SDup(i,j) = SPwavefunction(r(i,0),r(i,1), 2*j);
            SDdown(i,j) = SPwavefunction(r(i+nOrbitals,0),r(i+nOrbitals,1), 2*j+1);
        }
    }
    mat SDupinv = inv(SDup);
    mat SDdowninv = inv(SDdown);


    //Dot product of gradients
    mat JastrowGrad = zeros<mat>(nElectrons,nDims);
    mat SlaterGrad = zeros<mat>(nElectrons,nDims);


    //Slater laplacian and Jastrow laplacian
    double slaterLaplace = SlaterLaplace(nOrbitals,r);
    double jastrowLaplacian = JastrowLaplacian(r);//0;


    // compute the jastrow gradient and slater gradient
    for (int i = 0; i < nElectrons; i++) {
        for (int k = 0; k < nDims; k++) {

            double Jsum = 0.0;
            double Ssum = 0.0;
            //double J2sum = 0.0;

            for (int j = 0; j < i; j++) {
                double rji = RelativeDistance(r,j,i);
                Jsum += ((r(i,k)-r(j,k))/(rji))*a/((1.0+beta*rji)*(1.0+beta*rji));
            }

            for (int j = i+1; j < nElectrons; j++) {
                double rij = RelativeDistance(r,i,j);
                Jsum -= ((r(j,k)-r(i,k))/(rij))*a/((1.0+beta*rij)*(1.0+beta*rij));
            }
//            for (int j = i+1; j < nElectrons; j++) {
//                double rij = RelativeDistance(r,i,j);
//                Jsum += (r(i,k) - r(j,k))*a/(rij*pow((1.0+beta*rij),2));
//                J2sum += a/(rij*pow((1.0+beta*rij),2));
//                J2sum -= (r(i,k)-r(j,k))*(r(i,k)-r(j,k))*((1/rij)+4*beta+3*beta*beta*rij)/(rij*rij*pow((1+beta*rij),4));
//            }

            JastrowGrad(i,k) = Jsum;
            //jastrowLaplacian += J2sum;



            if (i < nOrbitals)
            {
                for (int j = 0; j < nOrbitals; j++)
                {
                    Ssum += SPpsiDer(2*j,r(i,0),r(i,1),k)*SDupinv(j,i);
                }
            }

            else
            {
                for (int j = 0; j < nOrbitals; j++)
                {
                    Ssum += SPpsiDer(2*j+1,r(i,0),r(i,1),k)*SDdowninv(j,i-nOrbitals);
                }
            }

            SlaterGrad(i,k) = Ssum;
        }
    }


    double GradientProduct = 2*dot(SlaterGrad,JastrowGrad);// + accu(JastrowGrad%JastrowGrad);



    //Total
    double EKin = -0.5*(jastrowLaplacian + slaterLaplace + GradientProduct);


    // Set up potential energy, external potential + eventual electron-electron repulsion
    //First: external potential
    double EPot = 0;

    for (int i = 0; i < nElectrons; i++) {
        double r2 = 0;
        for (int k = 0; k < nDims; k++) {
            r2 += r(i,k)*r(i,k);
        }
        EPot += omega*omega*0.5*r2;
    }

    double rij = 0;
    int count = 0;
    // Add the electron-electron repulsion
    for (int i = 0; i < nElectrons-1; i++) {
        for (int j = i+1; j < nElectrons; j++) {

            double r_12 = 0;
            for (int k = 0; k < nDims; k++) {
                  r_12 += (r(j,k)-r(i,k))*(r(j,k)-r(i,k));
            }
            if (Z) {
                EPot += (1.0/(sqrt (r_12)));
            }
            rij += sqrt(r_12);
            count++;

        }
    }

    double LocalE = EKin + EPot;
    //Z can be zero to turn of Coulomb repulsion
    rij = rij/count;
    vec Evec = zeros<vec>(4);

    Evec(0) = LocalE;
    Evec(1) = EKin;
    Evec(2) = EPot;
    Evec(3) = rij;

    return Evec;
}

double MonteCarlo::JastrowLaplacian(mat r)
{
    //calculates laplacian of jastrow factor

    double Laplace = 0;


    for (int k = 0; k < nElectrons; k++) {
        double sum1 = 0;
        for (int j = 0; j < nElectrons; j++) {
            if (j != k) {
                double rkj = RelativeDistance(r,k,j);
                for (int i = 0; i < nDims; i ++) {
                    sum1 += (2/rkj)*(r(k,i)-r(j,i))*a/(rkj*(1+beta*rkj)*(1+beta*rkj));
                    sum1 += a/(rkj*(1+beta*rkj)*(1+beta*rkj));
                    sum1 -= (r(k,i)-r(j,i))*(r(k,i)-r(j,i))*a*((1/rkj)+4*beta+3*beta*beta*rkj)/(rkj*rkj*pow((1+beta*rkj),4));
                }
            }

        }
        Laplace += sum1;
    }



    for (int k = 0; k < nElectrons; k++) {
        double sum2 = 0;
        vec rk(nDims,fill::zeros);
        rk(0) = r(k,0);
        rk(1) = r(k,1);

        for (int j = 0; j < nElectrons; j++ ) {
            if (j != k ) {
                vec rj(nDims,fill::zeros);
                rj(0) = r(j,0);
                rj(1) = r(j,1);
                double rkj = RelativeDistance(r,k,j);
                for (int i = 0; i < nElectrons; i++) {
                    if ( i != k) {
                        vec ri(nDims,fill::zeros);
                        ri(0) = r(i,0);
                        ri(1) = r(i,1);
                        double rki = RelativeDistance(r,k,i);

                        double vecprod = dot(rk-ri,rk-rj)/(rki*rkj);

                        for (int m = 0; m < nDims; m++) {
                            double derprod = ((r(k,m)-r(j,m))*a/(rkj*(1+beta*rkj)*(1+beta*rkj)))*((r(k,m)-r(i,m))*a/(rki*(1+beta*rki)*(1+beta*rki)));
                            sum2 += derprod;
                        }

                        sum2 *= vecprod;
                    }
                }


            }
        }
        Laplace += sum2;
    }
    return Laplace;

}

double MonteCarlo::computeLocalEnergy(const mat &r)
{
    //The function calculates the local energy using analytical expressions

    //Initialize
    int nOrbitals = nElectrons/2;
    mat distance = zeros<mat>(nElectrons, nElectrons);

    mat SDup(nOrbitals,nOrbitals,fill::zeros);
    mat SDdown(nOrbitals,nOrbitals,fill::zeros);

    //fill A with Hermite polynomial wavefunctions
    for (int i = 0; i < nOrbitals; i++) {
        for (int j = 0; j < nOrbitals; j++) {
            SDup(i,j) = SPwavefunction(r(i,0),r(i,1), 2*j);
            SDdown(i,j) = SPwavefunction(r(i+nOrbitals,0),r(i+nOrbitals,1), 2*j+1);
        }
    }
    mat SDupinv = inv(SDup);
    mat SDdowninv = inv(SDdown);

    // Set up interparticle distance
    for (int i = 0; i < nElectrons-1; i++) {
        for(int j = i+1; j < nElectrons; j++){
            distance(i,j) = RelativeDistance(r, i, j);
            distance(j,i) = distance(i,j);
        }
    }


    //Dot product of gradients
    mat JastrowGrad = zeros<mat>(nElectrons,nDims);
    mat SlaterGrad = zeros<mat>(nElectrons,nDims);



    // compute the jastrow gradient and slater gradient
    for (int i = 0; i < nElectrons; i++) {
        for (int k = 0; k < nDims; k++) {

            double Jsum = 0.0;
            double Ssum = 0.0;

            for (int j = 0; j < i; j++) {
                double rji = RelativeDistance(r,j,i);
                Jsum += ((r(i,k)-r(j,k))/(rji))*a/((1+beta*rji)*(1+beta*rji));
            }

            for (int j = i+1; j < nElectrons; j++) {
                double rij = RelativeDistance(r,i,j);
                Jsum -= ((r(j,k)-r(i,k))/(rij))*a/((1+beta*rij)*(1+beta*rij));
            }
            JastrowGrad(i,k) = Jsum;



            if (i < nOrbitals)
            {
                for (int j = 0; j < nOrbitals; j++)
                {
                    Ssum += SPpsiDer(2*j,r(i,0),r(i,1),k)*SDupinv(j,i);
                }
            }

            else
            {
                for (int j = 0; j < nOrbitals; j++)
                {
                    Ssum += SPpsiDer(2*j+1,r(i,0),r(i,1),k)*SDdowninv(j,i-nOrbitals);
                }
            }

            SlaterGrad(i,k) = Ssum;
        }
    }


    double GradientProduct = 2*dot(SlaterGrad,JastrowGrad);

    //Slater laplacian and Jastrow laplacian
    double slaterLaplace = SlaterLaplace(nOrbitals,r);
    double jastrowLaplacian = JastrowLaplacian(r);

    //Total
    double EKin = -0.5*(jastrowLaplacian + slaterLaplace + GradientProduct);


    // Set up potential energy, external potential + eventual electron-electron repulsion
    double EPot = 0.5*omega*omega*accu(r%r);

    // Add the electron-electron repulsion
    double ERep = 0.0;
    for (int i = 0; i < nElectrons-1; i++) {
        for (int j = i+1; j < nElectrons; j++) {
                ERep += 1.0/(distance(i,j));
        }
    }

    double LocalE = EKin + EPot + Z*ERep;

    return LocalE;
}

vec MonteCarlo::computeLocalEnergy_num(const mat &r)
{
    //The functio calculates the local enegy using numerical derivatives

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

    int count = 0;
    //Contribution from electron-electron potential (repulsion)
    for (int i = 0; i < nElectrons-1; i++) {
        for (int j = i+1; j < nElectrons; j++) {
            r1r2 = 0.0;
            for (int k = 0; k < nDims; k++) {
                r1r2 += (r(i,k) - r(j,k) ) * (r(i,k) - r(j,k) );
            }
            u += sqrt(r1r2);
            Eelel += 1/(sqrt(r1r2));
            count++;
        }
    }
    double rij = u/count;
    Eelel = Z * Eelel; //Can set Z = 0 to turn off Coulomb repulsion. Otherwise, it is 1

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

    double LocalE = EPot + EKin + Eelel;

    rij = rij/count;
    vec Evec = zeros<vec>(4);

    Evec(0) = LocalE;
    Evec(1) = EKin;
    Evec(2) = EPot+Eelel;
    Evec(3) = rij;

    return Evec;
}




