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

double CLa(double a,double CL1){
  //Lift force coefficient
  //function of angle of attack
  return CL1*a;//math.sin(a) //from Hummell (2003)
}

double CL(double a,double CL0,double CL1){
  //Total lift force coefficient.
  return CL0 + CLa(a,CL1);
}

double CDa(double a,double CD1){
  //Drag force coefficient
  //function of angle of attack
  a0 = pi/45.0;
  //According to Hummell (2003) this is the angle for which minimal drag is
  //achieved. This actually comes from data.
  return CD1*(a-a0)*(a-a0); //Hummell (2003)
}

double CD(double a,double CD0,double CD1){
  //Total drag for coefficient
  return CD0 + CDa(a,CD1);
}

double CRx(double wx,double CR1x){
  //Torque in the x-body direction (roll moment)
  //functoon of x-body angular velocity (roll angle)
  return -CR1x*wx; //Hummell (2003) long flight
}

double CRz(double wz,double CR1z){
  //Torque in the x-body direction (roll moment)
  //function of z-body angular velocity (spin)
  return -CR1z*wz;  //Hummell (2003)
}

double CR(double wx,double wz,double CR1x,double CR1z){
  //Total x-body torque
  return CRx(wx,CR1x) + CRz(wz,CR1z);
}

double CMa(double a,double CM1a){
  //Torque in the y-body direction (pitch moment)
  //function of the angle of attack
  return CM1a*a;//really not a good parameterization
}

double CMy(double wy,double CM1y){
  //Torque in the y-body direction (pitch moment)
  //function of the pitch angle angular velocity
  return -CM1y*wy; //Hummell (2003) long flight
}

double CM(double a,double wy,double CM0,double CM1a,double CM1y){
  //Total y-body torque
  return CM0 + CMa(a,CM1a) + CMy(wy,CM1y);
}

double CNz(double wz,double CN1z){
  //Torque in the z-body direction (spin down moment)
  //function of the yaw angle angular velocity (spin about the z-body axis)
  return -CN1z*wz; //Hummell (2003) long flight
}

double CN(double wz,double CN1z){
  //Total z-body torque. Just an alias for CNz made for consistent notation
  return CNz(wz,CN1z);
}
