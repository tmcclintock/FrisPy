import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def make_plots(positions, n_times):
	print "\nNow running animation routines"

	#Parse out all of the coordinates from the positions array
	t = positions[:,0]
	x = positions[:,1]
	y = positions[:,2]
	z = positions[:,3]
	vx = positions[:,4]
	vy = positions[:,5]
	vz = positions[:,6]
	phi = positions[:,7]
	theta = positions[:,8]
	phiDot = positions[:,9]
	thetaDot = positions[:,10]
	gammaDot = positions[:,11]
	gamma = positions[:,12]

	#Find the indices where the disc is still in the air
	in_air_indexes = np.where(z>0)
	print "Useful indicies: ",in_air_indexes[0][0]," through ",in_air_indexes[0][-1]

	#Plot Z vs X
	#plot_Z_vs_X(x,z,in_air_indexes)

	#Plot phi(t) and theat(t)
	#plot_angles(phi,theta,t,in_air_indexes)

	#End the plotting
	return

def plot_Z_vs_X(x,z,ind):
	print "\tPlotting Z vs X"
	fig = plt.figure()
	ax1 = plt.subplot(111)
	ax1.plot(x[ind],z[ind])
	ax1.set_ylim(0,max(z[ind]))
	plt.show()
	return

def plot_angles(phi,theta,t,ind):
	print "\tPlotting phi(t) and theta(t)"
	radians = np.pi/180. #used to convert to radians
	degrees = 180./np.pi #used to convert to degrees

	fig = plt.figure()
	plt.subplots_adjust(hspace=0.001)

	ax1 = plt.subplot(211)
	ax1.plot(t[ind],phi[ind]*degrees)
	ax1.set_xticklabels([])
	ax1.set_ylim(-180,180)

	ax2 = plt.subplot(212)#,sharex=ax1)
	ax2.plot(t[ind],theta[ind]*degrees)
	ax2.set_ylim(-90,90)
	fig.tight_layout()
	plt.show()
	return