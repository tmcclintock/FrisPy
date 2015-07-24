from setup import setup,config
import ctypes
import driver_interface

def FrisPy_Animation(config_filename):
	#Parse the configuration file
	#First check for necessary files for the animation
	conf = config.read_config_file(config_filename)
	if 'params_file_name' not in conf:
		raise Exception("Missing "+\
				"params_file_name in config_name")
	if 'initial_conditions_file_name' not in conf:
		raise Exception("Missing "+\
				"initial_conditions_file_name"+\
				" in config_name")

	#Read in the parameters (aka coefficients) and the initial conditions
	params_file_name = conf['params_file_name']
	initial_conditions_file_name = conf['initial_conditions_file_name']
	params,initial_conditions = setup.read_params_and_initial_conditions\
				    (params_file_name,
				     initial_conditions_file_name)
	
	#Pass the parameters and conditions to the actual animation routine
	positions = driver_interface.get_positions(initial_conditions,params)
	print positions, " the positions"