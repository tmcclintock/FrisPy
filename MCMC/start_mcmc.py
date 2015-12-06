import emcee
import numpy as np
import mcmc_likelihood as likelihood

#This function sets up emcee and acts as a driver for the MCMC search
#This returns an array of all of the samples and the likelihoods
def perform_mcmc(names,flight_data,initial_conditions,coeffs):
	print "\nPerforming an MCMC search using data for:\n\t",names

        #Set up the parameters list
        #NOTE: this design needs to be changed and have the initial
        #angles added as parameters
        #then the angles can be re-added to the initial_conditions list later
        mc_params = [initial_conditions,flight_data]

        #Define the number of dimensions and walkers and number of steps
        ndim = len(coeffs)
        nwalkers = 3*ndim
        nburnin = 1 #500
        nsteps = 3*nburnin

        #Define the initial positions of the walkers around the coeffs

        #Set up the sampler
        sampler = emcee.EnsembleSampler(nwalkers,ndim,likelihood.lnlike,args=mc_params)
        
        #Run the sampler on the burnin
        #pos,prob,state = sampler.run_mcmc(p0,nburnin)

        #Reset the chain
        sampler.reset()

        #Run the sampler for real
        #pos,prob,state = sampler.run_mcmc(pos,nsteps)

        #Return the chain and the probabilities
	#return sampler.chain,sampler.lnprobability
        return 0,0
