import numpy as np
from scipy import interpolate.InterpolatedUnivariateSpline as IUS
import driver_interface_mcmc

#This is a likelihood function that minimizes the percent difference
#between the observed trajectory and the model trajectory
def lnlike(positions,params):
    #Extract the parameters
    initial_positions,coeffs = params

    #Make the model
    #model_trajectory,n_times = driver_interface_mcmc.get_positions(initial_positions,coeffs)

    #Isolate the trajectory that matches the data (i.e. get rid of angles)

    #Create a spline for each of the coordinates as a function of the time

    #Compare the model at each time in the data to the trajectory at that time

    #Sum over the percent differences
    
    return 0
