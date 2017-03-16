#include "hermite.h"
#include <math.h>

#include <iostream>
using namespace std;


using namespace std;

hermite::hermite()
{

}


hermite::hermite(int n, double Z)
{
	H = calcHermite(n, Z); 
}

double hermite::calcHermite(int n, double Z)
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
		cout << "Wrong spatial quantum number given. It must be an integer between (or including) 0 and 9. Aborting."; 
		cout << endl; 
		exit(1); 
	}

	return H;

}
