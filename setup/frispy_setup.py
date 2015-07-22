import numpy as np

def read_params_and_data(params_file_name,data_file_name):
	#First read in the initial parameters
	#params_file_name = 'initial_params.txt'
	params_path = '../params/'+params_file_name
	initial_params = np.genfromtxt(params_path,names=True)

	print "\nStarting with initial parameters:"
	for i in xrange(0,len(initial_params.dtype.names)):
		print initial_params.dtype.names[i], initial_params[initial_params.dtype.names[i]]
	print '' #a blank line

	#Now read in the flight data
	#data_file_name = 'test_data.txt'
	data_path = '../data/'+data_file_name
	flight_data = np.genfromtxt(data_path,names=True)

	print '\nFlight data recieved from:'
	print '\t',data_path
	print '\nData recieved for:'
	print '\t',flight_data.dtype.names

	#Return the dictionary objects for both the initial parameters
	#and the flight data
	return[initial_params,flight_data]