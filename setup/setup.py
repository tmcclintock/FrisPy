import numpy as np


#Return a dictionary of the contents of the parameter file
def read_params(params_file_name):
	import numpy as np

	#Read in the initial parameters
	params_path = 'params/'+params_file_name
	params = np.genfromtxt(params_path,names=True)
	params_out = []

	print '\nParameters from:'
	print '\t',params_path
	print 'Using parameters:'
	for i in xrange(0,len(params.dtype.names)):
		print '\t',params.dtype.names[i], params[params.dtype.names[i]]
		params_out.append(params[params.dtype.names[i]])
	print '' #a blank line

	return params_out

#Return a dictionary of the contents of the initial conditions file
def read_initial_conditions(initial_conditions_file_name):
	import numpy as np

	#Reda in the initial conditions
	ic_path = 'data/'+initial_conditions_file_name
	initial_conditions = np.genfromtxt(ic_path,names=True)
	ic_out = []

	print '\nInitial conditions recieved from:'
	print '\t',ic_path
	print 'Initial conditions are:'
	for i in xrange(0,len(initial_conditions.dtype.names)):
		print '\t',initial_conditions.dtype.names[i],\
		      initial_conditions[initial_conditions.dtype.names[i]]
		ic_out.append(initial_conditions[initial_conditions.dtype.names[i]])

	#Return the initial conditions
	return ic_out

#Return a dictionary of the contents of the data file
def read_flight_data(flight_data_file_name):
	import numpy as np

	#Read in the flight data
	flight_data_path = 'data/'+flight_data_file_name
	flight_data = np.genfromtxt(flight_data_path,names=True)

	print '\nFlight data recieved from:'
	print '\t',flight_data_path
	print '\nData recieved for:'
	print '\t',flight_data.dtype.names

	#Return the dictionary objects for
	#the flight data
	return flight_data

#Return the dictionaries of the contents of the params and data files.
def read_params_and_flight_data(params_file_name,flight_data_file_name):
	#Read the params and flight data
	params = read_params(params_file_name)
	flight_data = read_flight_data(flight_data_file_name)

	#Return the dictionary objects for both the parameters
	#and the flight data
	return[params,flight_data]

#Return the dictionaries of the contents of the params and initial conditions
def read_params_and_initial_conditions(params_file_name,ic_file_name):
	#Read the params and initial conditions
	params = read_params(params_file_name)
	initial_conditions = read_initial_conditions(ic_file_name)

	#Return the dictionary objects for both the parameters
	#and the flight data
	return[params,initial_conditions]
