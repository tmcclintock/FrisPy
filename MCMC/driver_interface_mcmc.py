import ctypes
import numpy as np
#import mcmc_likelihood as likelihood

#This function interfaces with the driver.so shared library
#It gathers all of the times and positions to pass back to
#the mcmc routines
def get_positions(initial_positions,coeffs):
	return 0

#This function sets up emcee and acts as a driver for the MCMC search
#This returns an array of all of the samples and the likelihoods
def perform_mcmc(names,flight_data,params):
	print np.shape(names)
	print "Performing an MCMC search using data for:",names
	return 0,0 
