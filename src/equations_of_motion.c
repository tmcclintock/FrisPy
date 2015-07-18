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

void equations_of_motion(double*positions,double*derivs,
			    double t,
			    double*params){
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
  Tnc[0][0] = c_tht,  Tnc[0][1] = s_phi*s_tht, Tnc[0][2] = c_phi*s_tht;
  Tnc[1][0] = 0,      Tnc[1][1] = c_phi,       Tnc[1][2] = -s_phi;
  Tnc[2][0] = -s_tht, Tnc[2][1] = s_phi*c_tht, Tnc[2][2] = c_phi*c_tht;

  //Find the velocity components in the body frame
  double v_C_z = vx*Tnc[2][0] + vy*Tnc[2][1] + vz*Tnc[2][2];
  double v_C_x = 0;
  
}
