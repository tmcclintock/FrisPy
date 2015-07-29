import ctypes
import numpy as np

#This function interfaces with the driver.so shared library
#It gathers all of the times and positions to pass back to
#the animation routines
def get_positions(initial_positions,coeffs):
	print "Inside of get_positions"

	#Change these to copies numpy arrays so that
	#it is gauranteed that they are contiguous in memory
	#This is a REQUIREMENT for arrays in c and isn't
	#fulfilled in python
	initial_positions = np.array(initial_positions, dtype=np.double).copy()
	coeffs = np.array(coeffs, dtype=np.double).copy()


	#Animations shouldn't have specified flight times,
	#but since driver() takes one we specify a flight time
	#that is way longer than anything we would expect
	#a reasonable throw to take
	flight_time = 8 #seconds

	#Number of time steps
	n_times = 1000

	#Load the driver library and its funcitons
	#Then specify the arguement and return types
	driver_lib = ctypes.cdll.LoadLibrary('src/driver.so')
	driver = driver_lib.driver
	cleanup = driver_lib.cleanup
	driver.restype = ctypes.POINTER(ctypes.c_double)
	driver.argtypes = [ctypes.POINTER(ctypes.c_double),\
		ctypes.POINTER(ctypes.c_double),\
		ctypes.c_double,ctypes.c_int,\
		ctypes.POINTER(ctypes.c_double)]
	cleanup.restype = (None)
	cleanup.argtypes = [ctypes.POINTER(ctypes.c_double)]

	#Create the array that holds all of the positions and times
	all_positions = np.zeros(1000*13,dtype=np.double)
	print np.shape(all_positions),all_positions[0:10]

	#Create instances of pointers pointing to the important arrays
	ip_out = initial_positions.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	co_out = coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	ap_out = coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

	#Call the driver
	driver(ip_out,co_out,flight_time,n_times,ap_out)
	all_positions = np.array(ap_out[0:1000*13])
	print np.shape(all_positions),all_positions[0:10]

	#Call the cleanup
	#cleanup(all_positions)
	
	#Return the entire position array
	return 0#all_positions