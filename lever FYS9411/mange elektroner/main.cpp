#include <iostream>
#include "montecarlo.h"
#include <random>

#include "omp.h"
#include "mpi.h"

#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    cout <<  "Starting Monte Carlo program. " << endl;

    // Initialize parallellization
    int NumberProcesses = 4;
    int MyRank, NumberMCsamples;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &NumberProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &MyRank);
    double StartTime = MPI_Wtime();

    //Initialize timer
    clock_t start, stop;

    //Set important MC parameters
    int nDims = 2;
    int nEl = 2;

    cout << "Number of electrons: " << nEl << endl;

    int nCycles = 1e7; //for VMC
    double omega = 1.0; //harm osc frequency
    double timestep = 0.01; //timestep for importance sampling
    int MCcycles = 10000; //for steepest descent

    //With Coulomb?
    //bool Z = 0; //No
    bool Z = 1; //Yes

    //With Jastrow?
    //double a = 0; //No
    double a = 1.0; //Yes

    //input for regular VMC and guess for steepest descent
    double alpha = 0.99;
    double beta = 0.4;

    //Initialize MC object
    MonteCarlo MC4(nEl,Z,a,MCcycles,omega,timestep); //for steepest descent
    MonteCarlo MC5(nEl,Z,a,nCycles,omega,timestep,alpha,beta); //for regular VMC

    //Make matrices and vectors for stepest descent
    vec params = zeros<vec>(nDims); //params will contain final alpha and beta after SD
    vec guess = zeros<vec>(nDims); //guess contains guess for alpha and beta

    //Make initial guess for alpha, beta for steepest descent
    guess(0) = alpha;
    guess(1) = beta;

    double gamma = 0.01; //step size for SD

    //Run cycles! not parallel
    start = clock();
    MC5.runMonteCarloInt(); //VMC, with imp sampling
    //params = MC4.steepestDescent(guess, gamma); //steepest descent!
    stop = clock();

    cout << "Time used: " << (stop-start)/CLOCKS_PER_SEC << "s. " << endl;

    //From SD:
    cout << "Guess for alpha: " << guess(0) << endl;
    cout << "Estimated alpha: " << params(0) << endl;
    cout << "Guess for beta: " << guess(1) << endl;
    cout << "Estimated beta: " << params(1) << endl;


    //Run VMC calculation in paralell
//    NumberMCsamples = nCycles;
//    MPI_Bcast (&NumberMCsamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    ofstream outfile;
//    outfile.open("../pc/Dokumenter/MATLAB/FYS9411/project2/energies_6el.txt", ios::out | ios::binary);


//    double TotalEnergy, TotalEnergySquared, LocalProcessEnergy, LocalProcessEnergy2;
//    LocalProcessEnergy = LocalProcessEnergy2 = 0.0;
//    MC5.MonteCarloSampling(NumberMCsamples, LocalProcessEnergy, LocalProcessEnergy2);
//    //  Collect data in total averages
//    MPI_Reduce(&LocalProcessEnergy, &TotalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&LocalProcessEnergy2, &TotalEnergySquared, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    // Print out results  in case of Master node, set to MyRank = 0
//    if ( MyRank == 0) {
//        double Energy = TotalEnergy/( (double)NumberProcesses);
//        outfile << setiosflags(ios::showpoint | ios::uppercase);
//        outfile << setw(15) << setprecision(8) << Energy << endl;
//    }
//    double EndTime = MPI_Wtime();
//    double TotalTime = EndTime-StartTime;
//    if ( MyRank == 0 )
//        cout << "Time = " <<  TotalTime  << " on number of processors: "  << NumberProcesses  << endl;
//    if (MyRank == 0)
//        outfile.close();
//    // End MPI
//    MPI_Finalize();

    return 0;
}
