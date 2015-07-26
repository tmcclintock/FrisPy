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
	initia_positions = np.array(initial_positions, dtype=np.double).copy()
	coeffs = np.array(coeffs, dtype=np.double).copy()

	#Animations shouldn't have specified flight times,
	#but since driver() takes one we specify a flight time
	#that is way longer than anything we would expect
	#a reasonable throw to take
	flight_time = 8 #seconds

	#Load the driver library and its funcitons
	#Then specify the arguement and return types
	driver_lib = ctypes.cdll.LoadLibrary('src/driver.so')
	driver = driver_lib.driver
	cleanup = driver_lib.cleanup
	driver.restype = ctypes.POINTER(ctypes.c_double)
	driver.argtypes = [ctypes.POINTER(ctypes.c_double),\
		ctypes.POINTER(ctypes.c_double),\
		ctypes.c_double]
	cleanup.restype = (None)
	cleanup.argtypes = [ctypes.POINTER(ctypes.c_double)]

	#Call the driver
	#all_positions = driver(ctypes.cast(initial_positions,ctypes.POINTER(ctypes.c_double)),\
	#	ctypes.cast(coeffs,ctypes.POINTER(ctypes.c_double)),\
#		flight_time)

	#Call the cleanup
	#cleanup(all_positions)
	
	#Return the entire position array
	return all_positions