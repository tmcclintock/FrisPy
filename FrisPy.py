#Tom McClintock
#This is the main program for the Frisbee simulator, FrisPy.
#It is simply a small script that combines the rest of
#the parts of the program.

#Import the Animation and MCMC routines
from Animation import FrisPy_Animation as Animation
from MCMC import FrisPy_MCMC as MCMC

#Obtain the command line arguments
import sys
argv = sys.argv

#Check if no options have been given
if len(argv)>3 or len(argv)<3:
	message = "Usage error: incorrect number of arguments. Please run with 'python FrisPy [mode] [config_name]'."
	raise Exception(message)
#sys.exit(message)

#Choose either Animation or MCMC, and run in that mode or else throw an error.
mode = argv[1]
config_filename = argv[2]
if mode.lower() == 'animation':
	print "Running in Animation mode."
	Animation.FrisPy_Animation(config_filename)
elif mode.lower() == 'mcmc':
	print "Running in MCMC mode."
	MCMC.FrisPy_MCMC(config_filename)
else:
	message = "Usage Error: incorrect mode. Please specify either 'Animation' or 'MCMC'."
	raise Exception(message)
