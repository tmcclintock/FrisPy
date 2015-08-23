import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import animation_routines

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
	print "\tUseful indicies: ",in_air_indexes[0][0]," through ",in_air_indexes[0][-1]

	#Plot Z vs X
	animation_routines.plot_Z_vs_X(x,z,in_air_indexes)

	#Plot X vs Y
	animation_routines.plot_X_vs_Y(x,y,in_air_indexes)

	#Plot phi(t) and theat(t)
	animation_routines.plot_angles(phi,theta,t,in_air_indexes)

	#Do a 3D plot
	animation_routines.plot_3D_trajectory(x,y,z,len(in_air_indexes[0]))

	#End the plotting
	return