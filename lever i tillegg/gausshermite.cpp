// OpenMP code for numerical integration with Gauss-Hermite quadrature
// Here, functions are defined to compute the Coulomb integral in a Hartree-Fock calculation 
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <omp.h>

#include "gausshermite.h"

using namespace std;

gaussHermite::gaussHermite()
{
	//Default: harmonic oscillator frequency is hbar omega = 1 a. u. 
    omega = 1;
}

gaussHermite::gaussHermite(double w)
{
    omega = w;
}

int gaussHermite::factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// Function to define the integrand! Calls function defining hermite polynomials, as well as the harmonic oscillator mapping functions 
double gaussHermite::Integrand(double x1, double y1, double x2, double y2, int p, int q, int r, int s)
{

  double A = pow((omega/pi),0.25);

  double Ap = A*A*(1/sqrt(pow(2,harm.map2nx(p))*factorial(harm.map2nx(p))))*(1/sqrt(pow(2,harm.map2ny(p))*factorial(harm.map2ny(p))));
  double Aq = A*A*(1/sqrt(pow(2,harm.map2nx(q))*factorial(harm.map2nx(q))))*(1/sqrt(pow(2,harm.map2ny(q))*factorial(harm.map2ny(q))));
  double Ar = A*A*(1/sqrt(pow(2,harm.map2nx(r))*factorial(harm.map2nx(r))))*(1/sqrt(pow(2,harm.map2ny(r))*factorial(harm.map2ny(r))));
  double As = A*A*(1/sqrt(pow(2,harm.map2nx(s))*factorial(harm.map2nx(s))))*(1/sqrt(pow(2,harm.map2ny(s))*factorial(harm.map2ny(s))));

  double Hp = H.calcHermite(harm.map2nx(p),x1)*H.calcHermite(harm.map2ny(p),y1);
  double Hq = H.calcHermite(harm.map2nx(q),x2)*H.calcHermite(harm.map2ny(q),y2);;
  double Hr = H.calcHermite(harm.map2nx(r),x1)*H.calcHermite(harm.map2ny(r),y1);;
  double Hs = H.calcHermite(harm.map2nx(s),x2)*H.calcHermite(harm.map2ny(s),y2);;

  double z1 = x1/sqrt(omega);
  double z2 = x2/sqrt(omega);
  double w1 = y1/sqrt(omega);
  double w2 = y2/sqrt(omega);


  if  ((z1-z2)*(z1-z2)+(w1-w2)*(w1-w2) != 0)
    return Ap*Aq*Ar*As*Hp*Hq*Hr*Hs*1.0/(sqrt((z1-z2)*(z1-z2)+(w1-w2)*(w1-w2)));
  else
    return 0.0;
}

//Plain Gauss-Hermite integration with cartesian variables, brute force
double  gaussHermite::GaussHermiteIntegration(int n, int p, int q, int r, int s)
{
  double *x = new double [n];
  double *w = new double [n];
  GaussHermiteQuadrature(x, w, n);
  double Integral = 0.0;
  int i, j, k, l;
  # pragma omp parallel default(shared) private (i, j, k, l) reduction(+:Integral)
  {
    # pragma omp for
    for (i = 0;  i < n; i++){
      for (j = 0;  j < n; j++){
        for (k = 0;  k < n; k++){
            for (l = 0;  l < n; l++){
                Integral += w[i]*w[j]*w[k]*w[l]*Integrand(x[i],x[j],x[k],x[l],p,q,r,s);
            }
        }
      }
    }
  } // end parallel region
  delete [] x;
  delete [] w;
  return Integral;
}


// Setting Gaussian quadrature weights and integration points
void gaussHermite::GaussHermiteQuadrature(double *x, double *w, int n)
{
  int i,its,j,m;
  double p1,p2,p3,pp,z,z1;
  double Epsilon = 3.0e-14, PIM4 = 0.7511255444649425;
  int MaxIterations = 10;
  m=(n+1)/2;
  for (i=1;i<=m;i++) {
    if (i == 1) {
      z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
    } else if (i == 2) {
      z -= 1.14*pow((double)n,0.426)/z;
    } else if (i == 3) {
      z=1.86*z-0.86*x[0];
    } else if (i == 4) {
      z=1.91*z-0.91*x[1];
    } else {
      z=2.0*z-x[i-3];
    }
    for (its=1;its<=MaxIterations;its++) {
      p1=PIM4;
      p2=0.0;
      for (j=1;j<=n;j++) {
    p3=p2;
    p2=p1;
    p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
      }
      pp=sqrt((double)2*n)*p2;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= Epsilon) break;
    }
    if (its > MaxIterations) cout << "too many iterations in Hermite quadrature" << endl;
    x[i-1]=z;
    x[n-i] = -z;
    w[i-1]=2.0/(pp*pp);
    w[n-i]=w[i-1];
  }
} // Gaussian quadrature weights and integration points





void gaussHermite::GaussLegendreQuadrature(double x1, double x2, double *x, double *w, int n)
{
  int         m,j,i;
  double      z1,z,xm,xl,pp,p3,p2,p1;
  double      const  pi = 3.14159265359;
  double      *x_low, *x_high, *w_low, *w_high;
  double zero = 1.0e-10;

  m  = (n + 1)/2;                             // roots are symmetric in the interval
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  x_low  = x;                                       // pointer initialization
  x_high = x + n - 1;
  w_low  = w;
  w_high = w + n - 1;

  for(i = 1; i <= m; i++) {                             // loops over desired roots
    z = cos(pi * (i - 0.25)/(n + 0.5));

    /*
    ** Starting with the above approximation to the ith root
    ** we enter the mani loop of refinement bt Newtons method.
    */

    do {
      p1 =1.0;
      p2 =0.0;

      /*
      ** loop up recurrence relation to get the
      ** Legendre polynomial evaluated at x
      */

      for(j = 1; j <= n; j++) {
    p3 = p2;
    p2 = p1;
    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
      }

      /*
      ** p1 is now the desired Legrendre polynomial. Next compute
      ** ppp its derivative by standard relation involving also p2,
      ** polynomial of one lower order.
      */

      pp = n * (z * p1 - p2)/(z * z - 1.0);
      z1 = z;
      z  = z1 - p1/pp;                   // Newton's method
    } while(fabs(z - z1) > zero);

    /*
    ** Scale the root to the desired interval and put in its symmetric
    ** counterpart. Compute the weight and its symmetric counterpart
    */

    *(x_low++)  = xm - xl * z;
    *(x_high--) = xm + xl * z;
    *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
    *(w_high--) = *(w_low++);
  }
} // End_ function GaussLegendreQuadrature()




void gaussHermite::GaussLaguerreQuadrature(double *x, double *w, int n, double alf){

  int i,its,j;
  double ai;
  double p1,p2,p3,pp,z,z1;
  double Epsilon = 3.0e-14;
  int MaxIterations = 10;

  for (i=1;i<=n;i++) {
    if (i == 1) {
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {
      ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
        (1.0+3.5*ai))*(z-x[i-3])/(1.0+0.3*alf);
    }
    for (its=1;its<=MaxIterations;its++) {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
    p3=p2;
    p2=p1;
    p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      pp=(n*p1-(n+alf)*p2)/z;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= Epsilon) break;
    }
    if (its > MaxIterations) cout << "Too many iterations in Gauss Laguerre Quad" << endl;
    x[i-1]=z;
    w[i-1] = -exp(GammaFunction(alf+n)-GammaFunction((double)n))/(pp*n*p2);
  }
}


double gaussHermite::GammaFunction( double xx){
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


//  this function defines the function to integrate
double gaussHermite::int_function(double x)
{
  double value = x*x*x*exp(-x);
  return value;
} // end of function to evaluate



