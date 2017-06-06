#include <iostream>
#include "montecarlo.h"
#include <random>

#include <cmath>
#include <fstream>
#include <iomanip>


using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World! This is the program to run the MC cycles. " << endl;

    //Initialize timer
    clock_t start, stop;

    //Set important MC parameters
    int nDims = 2; //dimensions, must be 2
    int nEl = 2; //electrons, must be 2
    double steplength = 1.5; //for brute force VMC
    int nCycles = 5e5; //number of MC cycles
    double omega = 1.0; //frequency
    double timestep = 0.01; //for importance sampling
    int MCcycles = 10000; //for steepest descent

    //With Coulomb?
    //int Z = 0.0; //No
    int Z = 1.0; //Yes

    //With Jastrow?
    //double a = 0.0; //No
    double a = 1.0; //Yes

    //Set variational parameters
    double alpha = 0.990281;
    double beta = 0.399237;

    //Initialize MC object
    MonteCarlo MC; //with standard conditions
    MonteCarlo MC2(nEl, steplength,Z,a,nCycles,omega); //for the brute force version
    MonteCarlo MC3(nEl,Z,a,nCycles,omega,timestep); //for importance sampling
    MonteCarlo MC4(nEl,Z,a,MCcycles,omega,timestep); //for steepest descent
    MonteCarlo MC5(nEl,Z,a,nCycles,omega,timestep,alpha,beta);

    //Make matrices and vectors for stepest descent
    vec params = zeros<vec>(nDims); //params will contain alpha and beta
    vec guess = zeros<vec>(nDims); //guess contains guess for alpha and beta

    guess(0) = alpha;  //alpha
    guess(1) = beta ; //beta
    double gamma = 0.0001; //step for steepest descent

    //Run cycles!
    start = clock();

    //MC3.runMonteCarloInt(); //with imp sampling
    //MC2.runMonteCarloInt_bruteForce(); //brute force
    //params = MC4.steepestDescent(guess, gamma); //steepest descent!
    MC5.runMonteCarloInt();
    stop = clock();

    cout << "Time used: " << (stop-start)/CLOCKS_PER_SEC << "s. " << endl;

    cout << "Guess for alpha: " << guess(0) << endl;
    cout << "Estimated alpha: " << params(0) << endl;
    cout << "Guess for beta: " << guess(1) << endl;
    cout << "Estimated beta: " << params(1) << endl;

    return 0;
}
