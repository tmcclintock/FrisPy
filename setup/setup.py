import numpy as np

def read_params(params_file_name):
	#Read in the initial parameters
	#params_file_name = 'initial_params.txt'
	params_path = '../params/'+params_file_name
	initial_params = np.genfromtxt(params_path,names=True)

	print '\nStarting with initial parameters:'
	for i in xrange(0,len(initial_params.dtype.names)):
		print initial_params.dtype.names[i], initial_params[initial_params.dtype.names[i]]
	print '' #a blank line

	return initial_params

def read_data(data_file_name):
	#Read in the flight data
	data_path = '../data/'+data_file_name
	flight_data = np.genfromtxt(data_path,names=True)

	print '\nFlight data recieved from:'
	print '\t',data_path
	print '\nData recieved for:'
	print '\t',flight_data.dtype.names

	#Return the dictionary objects for both the initial parameters
	#and the flight data
	return flight_data

def read_params_and_data(params_file_name,data_file_name):
	#params_file_name = 'initial_params.txt'
	#data_file_name = 'test_data.txt'

	#Read the initial params and flight data
	initial_params = read_params(params_file_name)
	flight_data = read_data(data_file_name)

	#Return the dictionary objects for both the initial parameters
	#and the flight data
	return[initial_params,flight_data]