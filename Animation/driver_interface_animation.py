import ctypes
import numpy as np

#This function interfaces with the driver.so shared library
#It gathers all of the times and positions to pass back to
#the animation routines
def get_positions(initial_positions,coeffs):
	print "Inside of get_positions"

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
	#all_positions = driver(initial_positions,coeffs,flight_time)

	#Call the cleanup
	#cleanup(all_positions)
	
	#Return the entire position array
	return 0