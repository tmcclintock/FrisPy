from setup import setup,config
import ctypes
import numpy as np
import start_mcmc
#import mcmc_analysis as analysis

def FrisPy_MCMC(config_filename):
	#Parse config file
	conf = config.read_config_file(config_filename)
	if 'params_file_name' not in conf:
		raise Exception("Missing "+\
				"params_file_name in config_name")
	if 'flight_data_file_name' not in conf:
		raise Exception("Missing "+\
				"flight_data_file_name in config_name")

	#Read in the flight data
	params_file_name = conf['params_file_name']
	flight_data_file_name = conf['flight_data_file_name']
	params,flight_data = setup.read_params_and_flight_data(\
		params_file_name,flight_data_file_name)
	
	#flight_data is a numpy.ndarray with names,
        #while params is an array, so the flight
	#data names need to be extracted
	names = flight_data.dtype.names

	#Perform the MCMC
	samples,likelihoods = start_mcmc.perform_mcmc(\
		names,flight_data,params)

	#Create plots
	#analyze_samples.analysis(samples,likelihoods,names,flight_data,params)

	print "\nCompleting MCMC\n"
	return
