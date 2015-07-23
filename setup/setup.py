#Return a dictionary of the contents of the parameter file
def read_params(params_file_name):
	import numpy as np

	#Read in the initial parameters
	params_path = 'params/'+params_file_name
	params = np.genfromtxt(params_path,names=True)

	print '\nUsing parameters:'
	for i in xrange(0,len(params.dtype.names)):
		print params.dtype.names[i], params[params.dtype.names[i]]
	print '' #a blank line

	return params

#Return a dictionary of the contents of the initial conditions file
def read_initial_conditions(initial_conditions_file_name):
	import numpy as np

	#Reda in the initial conditions
	ic_path = 'data/'+initial_conditions_file_name
	initial_conditions = np.genfromtxt(ic_path,names=True)

	print '\nInitial conditions recieved from:'
	print '\t',ic_path
	print '\nInitial conditions for:'
	print '\t',initial_conditions.dtype.names

	#Return the initial conditions
	return initial_conditions

#Return a dictionary of the contents of the data file
def read_data(data_file_name):
	import numpy as np

	#Read in the flight data
	data_path = 'data/'+data_file_name
	flight_data = np.genfromtxt(data_path,names=True)

	print '\nFlight data recieved from:'
	print '\t',data_path
	print '\nData recieved for:'
	print '\t',flight_data.dtype.names

	#Return the dictionary objects for
	#the flight data
	return flight_data

#Return the dictionaries of the contents of the params and data files.
def read_params_and_data(params_file_name,data_file_name):
	#Read the params and flight data
	params = read_params(params_file_name)
	flight_data = read_data(data_file_name)

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
