import numpy as np
from scipy import interpolate.InterpolatedUnivariateSpline as IUS

#This is a likelihood function that minimizes the percent difference
#between the observed trajectory and the model trajectory
def lnlike():
    #Extract the parameters

    #Make the model

    #Isolate the trajectory that matches the data (i.e. get rid of angles)

    #Create a spline for each of the coordinates as a function of the time

    #Compare the model at each time in the data to the trajectory at that time

    #Sum over the percent differences
    
    return 0
