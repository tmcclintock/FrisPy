import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

def plot_Z_vs_X(x,z,ind):
	print "\tPlotting Z vs X"
	fig = plt.figure()
	ax1 = plt.subplot(111)
	ax1.plot(x[ind],z[ind])
	ax1.set_ylim(0,max(z[ind]))
	ax1.set_xlabel("X (m)")
	ax1.set_ylabel("Z (m)")
	plt.show()
	return

def plot_X_vs_Y(x,y,ind):
	print "\tPlotting Y vs X"
	fig = plt.figure()
	ax1 = plt.subplot(111)
	ax1.plot(x[ind],y[ind])
	#ax1.set_ylim(0,max(np.fabs(y[ind])))
	ax1.set_xlabel("X (m)")
	ax1.set_ylabel("Y (m)")
	plt.show()
	return

def plot_angles(phi,theta,t,ind):
	print "\tPlotting phi(t) and theta(t)"
	max_phi = max(np.fabs(phi[ind]))
	max_theta = max(np.fabs(theta[ind]))

	fig = plt.figure()
	plt.subplots_adjust(hspace=0.001)

	ax1 = plt.subplot(211)
	ax1.plot(t[ind],phi[ind])
	ax1.set_xticklabels([])
	ax1.set_ylim(-max_phi,max_phi)
	ax1.set_xlabel("Time (s)")
	ax1.set_ylabel(r"$\phi(t)$")

	ax2 = plt.subplot(212)#,sharex=ax1)
	ax2.plot(t[ind],theta[ind])
	ax2.set_ylim(-max_theta,max_theta)
	ax2.set_xlabel("Time (s)")
	ax2.set_ylabel(r"$\theta(t)$")

	fig.tight_layout()
	plt.show()
	return

def update_3D_trajectory(num,traj_lines,lines):
	for line, data in zip(lines,traj_lines):
		line.set_data(data[0:2,:num])
		line.set_3d_properties(data[2,:num])
	return line

#This creates a 3D trajectory of the throw
def plot_3D_trajectory(x,y,z,nsteps):
	print "\tPlotting 3D trajectory"
	mpl.rcParams['legend.fontsize'] = 10
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.set_xlim3d(min(x),max(x))
	ax.set_ylim3d(min(y),max(y))
	ax.set_zlim3d(0,max(z))
	fig.tight_layout()
	#The next line plots the trajectory immediately
	#full_line = ax.plot(x,y,z,label="Flight path",c='pink')
	lines = [ax.plot(x,y,z,label="Flight Path",c="b")[0],\
			 ax.plot(x,y,z*0,c="g",linestyle='--')[0],\
			 ax.plot(x,y*0+max(y),z,c="r",linestyle=":")[0]]
	ax.set_xlabel('X (m)')
	ax.set_ylabel('Y (m)')
	ax.set_zlabel('Z (m)')
	trajectory = np.array([[x,y,z],[x,y,z*0],[x,y*0+max(y),z]])
	lines_animation = animation.FuncAnimation(fig,update_3D_trajectory,\
			nsteps,fargs=(trajectory,lines),interval=0.0,blit=False)
	#lines_animation.save('example_throw.gif',fps=30)
	plt.show()
	return
