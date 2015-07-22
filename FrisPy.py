#Tom McClintock
#This is the main program for the Frisbee simulator, FrisPy.
#It is simply a small script that combines the rest of
#the parts of the program.

#Obtain the command line arguments
import sys
argv = sys.argv

#Check if no options have been given
if len(argv)>2 or len(argv)<2:
	message = "Usage error: too many options. Please Specify either Animation or MCMC."
	sys.exit(message)

mode = argv[1]
print mode

if mode is 'a':
	print mode