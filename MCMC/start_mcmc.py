import emcee
import numpy as np
import mcmc_likelihood as likelihood

#This function sets up emcee and acts as a driver for the MCMC search
#This returns an array of all of the samples and the likelihoods
def perform_mcmc(names,flight_data,initial_conditions,params):
	print np.shape(names)
	print "Performing an MCMC search using data for:",names

        mc_params = [initial_conditions,params]

	return 0,0 
