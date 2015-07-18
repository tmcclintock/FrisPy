#include <stdio.h>
#include <math.h>

#include "coefficients.h"
#include "constants.h"

/* Naming convention:
   - any variable with a 'Dot' at the end means it is a time derivative
   - 'positions' refers to x, y, z as well as the euler angles in addition
     to the time derivatives of all of them
   - 'derivs' are the time derivatives of all of the positions
   - 't' has the current time, which is not needed in this calculation,
     but could in principle account for a time-varying wind
   - 'params' contains all of the parameters we are fitting for. That is,
     all of the flight parameters that are unknown
*/

double dot(double*A,double*B){
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void cross(double*A,double*B,double*C){
  //Uses a temp array so that cross produces can be done in place
  double D[3];
  D[0] = A[1]*B[2]-A[2]*B[1];
  D[1] = A[2]*B[0]-A[0]*B[2];
  D[2] = A[0]*B[1]-A[1]*B[0];
  C[0]=D[0],C[1]=D[1],C[2]=D[2];
  return;
}

void add(double*A,double *B,double *C){
  //Uses a temp array so that cross produces can be done in place
  double D[3];
  D[0] = A[0]+B[0];
  D[1] = A[1]+B[1];
  D[2] = A[2]+B[2];
  C[0]=D[0],C[1]=D[1],C[2]=D[2];
  return;
}

void equations_of_motion(double*positions,double*derivs,
			    double t,
			    double*coeffs){
  //Declare some temporary variables that may be used
  double temp1, temp2;
  int i,j,k;

  /*The six coordinates and velocities in the lab frame
    arrive in the arrays 'positions'.
   */
  double x = positions[0], y = positions[1], z = positions[2];
  double vx= positions[3], vy= positions[4], vz= positions[5];
  double phi     = positions[6], theta    = positions[7];
  double phiDot  = positions[8], thetaDot = positions[9];
  double gammaDot= positions[10], gamma = positions[11];

  //Take in the various force and moment coefficients
  double CL0 = coeffs[0], CL1 = coeffs[1];
  double CD0 = coeffs[2], CD1 = coeffs[3];
  double CR1x = coeffs[4], CR1z = coeffs[5];
  double CM0 = coeffs[6], CM1a = coeffs[7], CM1y = coeffs[8];
  double CN1z = coeffs[9];

  //Check if the disc has hit the ground
  if (z <= 0){
    for(i=0; i<12; i++)
      derivs[i] = 0;
    return;
  }

  //TODO: include a wind. Probably as a static set of variables.
  /*
    vx = vx + v_wind_x;
    vy = vy + v_wind_y;
    vz = vz + v_wind_z;
  */

  //Calculate the sine and cosine of theta and phi to be used later
  //This is a convenient short hand
  double s_phi = sin(phi),  c_phi = cos(phi);
  double s_tht = sin(theta),c_tht = cos(theta);

  //Construct the Euler rotation matrix to go from 
  //the inertial frame (N) to the body frame (C)
  // the 
  double Tnc[3][3];
  Tnc[0][0] = c_tht,  Tnc[0][1] = s_phi*s_tht, Tnc[0][2] = -c_phi*s_tht;
  Tnc[1][0] = 0,      Tnc[1][1] = c_phi,       Tnc[1][2] = s_phi;
  Tnc[2][0] = s_tht,  Tnc[2][1] = -s_phi*c_tht, Tnc[2][2] = c_phi*c_tht;

  //Define the rotation matrix rows
  double C1[] = {Tnc[0][0],T[0][1],T[0][2]};
  double C2[] = {Tnc[1][0],T[1][1],T[1][2]};
  double C3[] = {Tnc[2][0],T[2][1],T[2][2]};

  //Find the velocity components in the body frame
  double v[]={vx,vy,vz};
  double v_dot_C3 = dot(v,C3);
  double v_plane[] = {v[0]*(1-v_dot_C3),v[1]*(1-v_dot_C3),v[2]*(1-v_dot_C3)};

  //Find the length of these variace velocities
  double norm_v_plane = sqrt(dot(v_plane,v_plane));
  double norm_v = sqrt(dot(v,v));
  double norm_vp = sqrt(dot(v_plane,v_plane));

  //Find the unit vectors for various directions
  double v_hat[] = {v[0]/norm_v,v[1]/norm_v,v[2]/norm_v};
  double vp_hat[] = {v_plane[0]/norm_vp,v_plane[1]/norm_vp,v_plane[2],norm_vp}
  double ulat[3];
  cross(C3,vp_hat,ulat);

  //Calculate the angle of attack, alpha
  //This angle is negative since if there is face on the disc 
  //it should have a positive angle of attack even though vz_C
  //would be negative.
  double alpha = -atan(v_dot_C3/norm_v_plane);

  //Make copies of some of the arrays, mostly so that it is easier to see what is happening
  double x_C_hat[] = {vp_hat[0],vp_hat[1],vp_hat[2]};
  double y_C_hat[] = {ulat[0],ulat[1],ul[2]};
  double z_C_hat[] = {C3[0],C3[1],C3[2]};

  //Calculate the lift force
  double Flift_amp = CL(alpha,CL0,CL1)*A*rho*dot(v,v)/2.;
  double Flift[3];
  cross(v_hat,y_C_hat,Flift);
  Flift[0] = Flift[0]*Flift_amp; Flift[1] = Flift[1]*Flift_amp; Flift[2] = Flift[2]*Flift_amp;

  //Calculate the drag force
  
}
