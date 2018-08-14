#include <math.h>
#include <stdio.h>

#include "equations_of_motion.h"

void rk4(double*positions,double dt,double t,double*coeffs){
	//Declare the K arrays that hold the derivatives
	double k1[12];
	double k2[12];
	double k3[12];
	double k4[12];

	//A temporary array
	double temp[12];

	//Iteration variables
	int i,j;

	//If the disk has hit the ground, then keep the same positions
	if(positions[2] <= 0){
		//if z <= 0
		return;
	}

	//Find the derivatives at the start time for k1
	equations_of_motion(positions,k1,t,coeffs);

	//Calculate the temp array
	for(i=0;i<12;i++){
		temp[i] = positions[i]+dt/2.*k1[i];
		printf("%e\t",k1[i]);
	}
	printf("\n");

	//Find the derivatives at t+dt/2 for k2
	equations_of_motion(temp,k2,t+dt/2.,coeffs);

	//Calculate the new temp arrays
	for(i=0;i<12;i++)
		temp[i] = positions[i]+dt/2.*k2[i];

	//Find the derivatives at t+dt/2 for k3
	equations_of_motion(temp,k3,t+dt/2.,coeffs);

	//Calculate the new temp arrays
	for(i=0;i<12;i++)
		temp[i] = positions[i]+dt*k3[i];

	//Find the derivatives at t+dt for k4
	equations_of_motion(temp,k4,t+dt,coeffs);

	//Update the positions
	for (i=0;i<12;i++){
		positions[i] = positions[i] + dt/6.*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i]);
	}
	//End rk4
	return;
}
