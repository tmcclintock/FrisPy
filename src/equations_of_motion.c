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
  -  'tau' refers to a torque (aka moment) along one of the inertial axes
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
  //NOTE: Tnc takes you from N to C
  //while {Tnc}^T = Tcn takes you from C to N
  double Tnc[3][3];
  Tnc[0][0] = c_tht, Tnc[0][1] = s_phi*s_tht,  Tnc[0][2] = -c_phi*s_tht;
  Tnc[1][0] = 0,     Tnc[1][1] = c_phi,        Tnc[1][2] = s_phi;
  Tnc[2][0] = s_tht, Tnc[2][1] = -s_phi*c_tht, Tnc[2][2] = c_phi*c_tht;

  //Define the rotation matrix rows
  double C1[] = {Tnc[0][0],Tnc[0][1],Tnc[0][2]};
  double C2[] = {Tnc[1][0],Tnc[1][1],Tnc[1][2]};
  double C3[] = {Tnc[2][0],Tnc[2][1],Tnc[2][2]};

  //Find the velocity components in the body frame
  double v[]={vx,vy,vz};
  double v_dot_C3 = dot(v,C3);
  double v_plane[] = {v[0]*(1-v_dot_C3),v[1]*(1-v_dot_C3),v[2]*(1-v_dot_C3)};

  //Find the length of these various velocities
  double norm_v_plane = sqrt(dot(v_plane,v_plane));
  double norm_v = sqrt(dot(v,v));
  double norm_vp = sqrt(dot(v_plane,v_plane));

  //Find the unit vectors for various directions
  double v_hat[] = {v[0]/norm_v,v[1]/norm_v,v[2]/norm_v};
  double vp_hat[] = {v_plane[0]/norm_vp,v_plane[1]/norm_vp,v_plane[2],norm_vp};
  double ulat[3];
  cross(C3,vp_hat,ulat);

  //Calculate the angle of attack, alpha
  //This angle is negative since if there is face on the disc 
  //it should have a positive angle of attack even though vz_C
  //would be negative.
  double alpha = -atan(v_dot_C3/norm_v_plane);

  //Make copies of some of the arrays, mostly so that it is easier to see what is happening
  double x_C_hat[] = {vp_hat[0],vp_hat[1],vp_hat[2]};
  double y_C_hat[] = {ulat[0],ulat[1],ulat[2]};
  double z_C_hat[] = {C3[0],C3[1],C3[2]};

  //Calculate the lift force
  double Flift_amp = CL(alpha,CL0,CL1)*Area*rho*dot(v,v)/2.;
  double Flift[3];
  cross(v_hat,y_C_hat,Flift);
  Flift[0] = Flift[0]*Flift_amp; Flift[1] = Flift[1]*Flift_amp; Flift[2] = Flift[2]*Flift_amp;

  //Calculate the drag force
  double Fdrag_amp = CD(alpha,CD0,CD1)*Area*rho*dot(v,v)/2.;
  double Fdrag[]={-Fdrag_amp*v_hat[0],-Fdrag_amp*v_hat[1],-Fdrag_amp*v_hat[2]};

  //Calculate the gravitational force
  double Fgrav[]={0,0,-m*g};

  //Sum the forces to get the total force
  double Ftot[] = {Flift[0]+Fdrag[0],Flift[1]+Fdrag[1],Flift[2]+Fdrag[2]+Fgrav[2]};

  //Write the angular velocities in the body frame
  double w_in_C[] = {phiDot*c_tht, thetaDot, phiDot*s_tht + gammaDot};

  //Calculate the angular velocity in the inertial frame
  double w_in_N[3];
  for(i=0;i<3;i++){
    temp1 = 0;
    for(j=0;j<3;j++){//{Tnc}^T dot w_in_C
      temp1+= Tnc[j][i] * w_in_C[j];
    }
    w_in_N[i] = temp1;
  }

  //Find the components of the angular velocity
  //as expressed in the inertial frame along each of the body frame unit
  //vectors
  double w_x = dot(w_in_N,x_C_hat);
  double w_y = dot(w_in_N,y_C_hat);
  double w_z = dot(w_in_N,z_C_hat);

  //Calculate the spin parameter, which could be used later
  //to incorporate the Robins-Magnus forec or if we
  //have more accurate data about drag coefficients
  double spin_parameter = w_z*d/(2*norm_v);

  //Calculate each torque (aka moment)
  double tau_amp = Area*rho*d*dot(v,v)/2;
  double tau_x_amp = CR(w_x,w_z,CR1x,CR1z)*tau_amp;
  double tau_y_amp = CM(alpha,w_y,CM0,CM1a,CM1y)*tau_amp;
  double tau_z_amp = CN(w_z,CN1z)*tau_amp;
  double tau_x[] = {x_C_hat[0]*tau_x_amp, x_C_hat[1]*tau_x_amp, x_C_hat[2]*tau_x_amp};
  double tau_y[] = {y_C_hat[0]*tau_y_amp, y_C_hat[1]*tau_y_amp, y_C_hat[2]*tau_y_amp};
  double tau_z_N[] = {0, 0, tau_z_amp};
  //Note: tau_x and tau_y are in the body frame, while tau_z is already in N

  //Calculate tau_x_N and tau_y_N
  double tau_x_N[3];
  double tau_y_N[3];
  for(i=0;i<3;i++){
    temp1 = 0;
    temp2 = 0;
    for(j=0;j<3;j++){//{Tnc}^T dot w_in_C
      temp1+= Tnc[i][j] * tau_x[j];
      temp2+= Tnc[i][j] * tau_y[j];
    }
    tau_x_N[i] = temp1;
    tau_y_N[i] = temp2;
  }

  //Calculate the total torque (aka moment) in the inertial frame
  double tau_total[] = {tau_x_N[0]+tau_y_N[0],
    tau_x_N[1]+tau_y_N[1],
    tau_x_N[2]+tau_y_N[2]+tau_z_N[2]};

  //If we want, set the torques all to 0
  //tau_total[0] = 0, tau_total[1] = 0, tau_total[2] = 0;

  //TODO: Remove the wind from the velocities
  /*
    vx = vx + v_wind_x;
    vy = vy + v_wind_y;
    vz = vz + v_wind_z;
  */

  //Calculate all of the derivatives
  derivs[0] = vx;
  derivs[1] = vy;
  derivs[2] = vz;
  derivs[3] = Ftot[0]/m;
  derivs[4] = Ftot[1]/m;
  derivs[5] = Ftot[2]/m;
  derivs[6] = phiDot;
  derivs[7] = thetaDot;
  derivs[8] = (tau_total[0] 
               + 2*Iss*thetaDot*phiDot*s_tht 
               - Izz*thetaDot*(phiDot*s_tht+gammaDot)
              )/(Iss*c_tht);
  derivs[9] = (tau_total[1]
               + Izz*phiDot*c_tht*(phiDot*s_tht+gammaDot)
               - Iss*phiDot*phiDot*c_tht*s_tht
              )/Iss;
  derivs[10]= (tau_total[2]
               - Izz*derivs[8]*s_tht
               - Iss*thetaDot*phiDot*c_tht
              )/Izz;
  derivs[11]= gammaDot;

  //End the function
  return;
}
