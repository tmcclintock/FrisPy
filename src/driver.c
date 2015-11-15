#include <stdio.h>
#include <stdlib.h>

#include "rk4.h"

/*
	initial_position: the initial configuration of the frisbee at the start of its flight
	coeffs: the coefficients used in the current calculation
	flight_time: the total time of flight
*/
void driver(void*initial_positionav,void*coeffsav,double flight_time,int n_times,double*all_positions){
	//The input arrays came in as void* so we need to cast them back
	double*initial_position = (double*)initial_positionav;
	double*coeffs = (double*)coeffsav;

	//Declare some iteration variables
	int i,j;

	//These are the the number of variables
	int n_vars=12;

	//Declare an array for the current positions
	double current_position[n_vars];

	//Calculate the time step
	double dt = flight_time/n_times;

	printf("\tTimestep = %e\n",dt);

	//Declare the current time
	double t=0;
	all_positions[0]=t;

	//Put in the current position as the initial position
	for (i = 0; i < n_vars; ++i)
	{
		current_position[i] = initial_position[i];
		all_positions[i+1] = current_position[i];
	}

	//Loop over all the times and advance the simulation
	//Note: we take n_times-1 steps since we have already stored the initial information
	for (i = 0; i < n_times-1; ++i)
	{
		rk4(current_position,dt,t,coeffs);
		t+=dt;

		//Store the current position and time
		all_positions[(i+1)*(n_vars+1)+0]=t;
		for (j = 0; j < n_vars; ++j)
		{
			all_positions[(i+1)*(n_vars+1)+(j+1)] = current_position[j];
		}
	}

	//Now return the all_positions array
	return;
}

//This function is necessary to free c arrays, since that can't be
//done from within python
void cleanup(double*array){
	free(array);
	return;
}
