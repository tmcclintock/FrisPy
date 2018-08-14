#include <stdio.h>
#include <math.h>

//Define coefficients that have functional forms
//These functional forms attempt to mimic the experimental
//data presented in Hummell (2003) gathered from studies such
//as Potts and Crowther (2002), Yasada (1999),
//and Stilley and Carstens (1972).
//The functional forms themselves are explained within the functions.
//Overall, they need to be periodic (obviously) and need to vary
//in sensible ways.

double CLa(double a,double CL1);
double CL(double a,double CL0,double CL1);
double CDa(double a,double CD1);
double CD(double a,double CD0,double CD1);
double CRx(double wx,double CR1x);
double CRz(double wz,double CR1z);
double CR(double wx,double wz,double CR1x,double CR1z);
double CMa(double a,double CM1a);
double CMy(double wy,double CM1y);
double CM(double a,double wy,double CM0,double CM1a,double CM1y);
double CNz(double wz,double CN1z);
double CN(double wz,double CN1z);
