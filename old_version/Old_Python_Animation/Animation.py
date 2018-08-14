#Written by Tom McClintock
#2015 University of Arizona Physics
#
#This program simulates the flight of frisbees **ahem** flying discs.
#It does this by considering drag, lift, and gravitational forces
#as well as inertial moments that affect the rotational rates of
#the disc about its three principle axes.
#
#If you cannot tell, this program is written in python. It uses
#the math, numpy, scipy, and matplotlib libraries which are all free to
#download.
#
#Computationally speaking this simulation is very accurate. It uses
#a 4th order runge-kutta integrator.
#
#By far, the most uncertainty comes from a lack of details of
#the lift, drag, and spin coefficients. Attempts have been made
#to measure some of these, and here we borrow largely from
#Hummell (2003). However, in order to make a more robust model
#we estimate a funcitonal form of these coefficients as a function
#of the body angles in such a way as to make them periodic from 0
#to 2pi, to reflect the rotational symmetry of the disc in both the
#pitch and the roll. The functional forms are ALMOST periodic
#from 0 to pi, to reflect the fact that a disc thrown upside down
#(e.g. a hammer) can exhibit sensible flight behaviour.
#
#Feel free to contact the author of the code at
#tmcclintock@physics.arizona.edu
#
#This code is very similar to Hummell (2003), and for that reason the
#author is greatly indebted to her for her work. However, it was
#decided that it would be worth rewriting her code
#in a language that wasn't proprietary. To Sarah Hummell, hopefully
#I meet you on the field on day!

####################################################################


#Libraries used:
#Math, Numpy, the ODE solvers in scipy, and pyplot in matplotlib
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import math
throwname = "fast"

####################################################################


#Define constants
pi = math.pi
m = 0.175         #kg; frisbee mass
g = 9.81          #m/s^2; acceleration due to gravity
A = 0.057         #m^2; planform (top down) area
d = 2*(A/pi)**0.5 #m; diameter
rho = 1.23        #kg/m^3; average air density
Iz = 0.002352     #kg*m^2; moment of inertia about zz
Ixy = 0.001219    #kg*m^2; moment of inertia about either xx or yy

####################################################################

#Declare a vector that holds the wind
vwind = np.array([0,0,0]) #m/s; three components of the wind velocity


####################################################################


#Define constant coefficients
#All of these numbers are from Hummell (2003)
CL0 = 0.3331   #Lift coefficient constant
CD0 = 0.1769   #Drag coefficient constant
CM0 = -0.0821  #Pitch angle torque (moment) coefficient constant

####################################################################

#Define our main funciton
def main():
    #Set the initial variable conditions
    #Long flat backhand
    speed0 = 14.0    #m/s; initial speed
    phi0  = 0.0   #degrees; roll angle (+ means OI BH/IO FH, - is opposite)
    tht0 = 0.0       #degrees; polar angle tilt (- means the throw is airbounce-y)
    fpAng0 = 0.0       #degrees; flight path angle between v and x-lab
    gamDot0 = -50.0    #rad/s; spin on the disc (- is BH, + is FH)
    xyang = 00.0     #angle between x-lab axis and velocity. This is useful for plotting
    tf = 3.2  #s; final time
    throwname = "fast2"


    #Gentle slightly rising forehand
    #speed0 = 9.0    #m/s; initial speed
    #phi0  = 0.0   #degrees; roll angle (+ means OI BH/IO FH, - is opposite)
    #tht0 = -10.0       #degrees; polar angle tilt (- means the throw is airbounce-y)
    #fpAng0 = 0.0       #degrees; flight path angle between v and x-lab
    #gamDot0 = -50.0    #rad/s; spin on the disc (- is BH, + is FH)
    #xyang = 0     #angle between x-lab axis and velocity. This is useful for plotting
    #tf = 2.2  #s; final time
    #throwname = "slow"

    #Similar to above, but should get obliterated by wind
    #speed0 = 9.0    #m/s; initial speed
    #phi0  = -10.0   #degrees; roll angle (+ means OI BH/IO FH, - is opposite)
    #tht0 = -10.0       #degrees; polar angle tilt (- means the throw is airbounce-y)
    #fpAng0 = 0.0       #degrees; flight path angle between v and x-lab
    #gamDot0 = -50.0    #rad/s; spin on the disc (- is BH, + is FH)
    #xyang = 0     #angle between x-lab axis and velocity. This is useful for plotting
    #tf = 3.2  #s; final time
    #throwname = "windy"
    #vwind[0] = 6 #x direction
    #vwind[1] = 0 #y direction
    #vwind[2] = 0 #z direction

    
    #Set the wind speed in each direction. Units are m/s
    vwind[0] = 0 #x direction
    vwind[1] = 0 #y direction
    vwind[2] = 0 #z direction

    #Note: a negative gamDot means the disc is spinning clockwise looking down
    #i.e. a backhand release. Thus a positive gamDot is a forehand.

    #Set the initial conditions to be fed to the integrator
    x = 0.  #m; x-lab displacement
    y = 0.  #m; y-lab displacement
    z = 1.  #m; height above the ground
    vx = speed0*math.cos(fpAng0*2*pi/360)*math.cos(xyang*2.0*pi/360.) #initial x-component velocity
    vy = speed0*math.sin(xyang*2.0*pi/360.)           #initial y-component velocity
    vz = speed0*math.sin(fpAng0*2*pi/360) #initial z-component velocity
    phi = phi0*2*pi/360#degrees; roll angle 
    tht = tht0*2*pi/360#degrees; pitch angle 
    phiDot = 0.0        #radians/s; roll angular momentum
    thtDot = 0.0        #radians/s; pitch anglular momentum
    gamDot = gamDot0          #radians/s; initial z-body angular momentum (spin)
    gam = 0.     #radians; spin angle (has no effect on computation)
    
    #Construct the initial parameter array
    params = np.array([x,y,z,vx,vy,vz,phi,tht,phiDot,thtDot,gamDot,gam])

    #Declare the initial and final times
    ti = 0.0   #s; initial time

    #Declare the number of time steps and find the step size
    nsteps = 100
    dt = (tf-ti)/nsteps

    #Declare the array that contains the times
    print("Timestep = %4.3f seconds" % dt)
    t = np.arange(ti,tf,dt)

    #Integrate the solution
    #solution is an nsteps X 12 size array containing the solution at all steps
    solution = odeint(eqOfMotion,params,t)

    #Print the solution
    #print(solution[:,0])

    #Plot the angles as a funciont of time
    plotAngles(solution,t,nsteps)

    #Plot the trajectory
    plotTrajectory(solution,nsteps)

####################################################################
####################################################################
####################################################################
####################################################################



#Define coefficients that have functional forms
#These functional forms attempt to mimic the experimental
#data presented in Hummell (2003) gathered from studies such
#as Potts and Crowther (2002), Yasada (1999),
#and Stilley and Carstens (1972).
#The functional forms themselves are explained within the functions.
#Overall, they need to be periodic (obviously) and need to vary
#in sensible ways.
def CLa(a):
    #Lift force coefficient
    #function of angle of attack
    return 1.9124*math.sin(a) #from Hummell (2003)

def CL(a):
    #Total lift force coefficient.
    return CL0 + CLa(a)


def CDa(a):
    #Drag force coefficient
    #function of angle of attack
    a0 = pi/45
    #According to Hummell (2003) this is the angle for which minimal drag is
    #achieved. This actually comes from data.
    return 0.685*(a-a0)*(a-a0) #Hummell (2003)

def CD(a):
    #Total drag for coefficient
    return CD0 + CDa(a)

def CRx(wx):
    #Torque in the x-body direction (roll moment)
    #functoon of x-body angular velocity (roll angle)
    return -0.0125*wx #Hummell (2003) long flight

def CRz(wz):
    #Torque in the x-body direction (roll moment)
    #function of z-body angular velocity (spin)
    return -0.00171*wz  #Hummell (2003)

def CR(wx,wz):
    #Total x-body torque
    return CRx(wx) + CRz(wz)

def CMa(a):
    #Torque in the y-body direction (pitch moment)
    #function of the angle of attack
    return 0.4338*math.sin(a)#really not a good parameterization

def CMy(wy):
    #Torque in the y-body direction (pitch moment)
    #function of the pitch angle angular velocity
    return -0.0144*wy #Hummell (2003) long flight

def CM(a,wy):
    #Total y-body torque
    return CM0 + CMa(a) + CMy(wy)

def CNz(wz):
    #Torque in the z-body direction (spin down moment)
    #function of the yaw angle angular velocity (spin about the z-body axis)
    return -0.0000341*wz #Hummell (2003) long flight

def CN(wz):
    #Total z-body torque. Just an alias for CNz made for consistent notation
    return CNz(wz)

####################################################################


#Define the differential equation for each variable
def eqOfMotion(params,t):
    #Acquire the input parameters
    #These are the variables for which we must calculate time derivatives for
    x,y,z,vx,vy,vz,phi,tht,phiDot,thtDot,gamDot,gam = params[0:12]

    #Declare the outgoing derivatives array. Each element aligns with the
    #incoming elements in params
    paramsDot = np.zeros(12)

    #If the disc has hit the ground, return all zeros for the derivatives
    if z<0:return paramsDot
    
    #Apply the wind to the velocities.
    #Note: this is undone later on, so as to not move the disc too far
    vx = vx + vwind[0]
    vy = vy + vwind[1]
    vz = vz + vwind[2]
    
    #Construct the Euler rotation matrix to go from the inertial frame (N)
    #to the body frame (C)
    Tnc = np.array(
        [[math.cos(tht),math.sin(phi)*math.sin(tht),-math.cos(phi)*math.sin(tht)],
        [0,math.cos(phi),math.sin(phi)],
        [math.sin(tht),-math.sin(phi)*math.cos(tht),math.cos(phi)*math.cos(tht)]])
    
    #These are the rows of Tnc.
    #C3 is also the z-body unit (hat) vector expressed in inertial coordinates
    C1 = Tnc[0]
    C2 = Tnc[1]
    C3 = Tnc[2]

    #Create the velocity vector in the lab frame
    v = np.array([vx,vy,vz])

    #Find v dot C3, then vp (this is the velocity in the plane of the disc)
    vC3 = np.dot(v,C3)
    vp = v - C3*vC3

    #Calculate the angle of attack, alpha (alf)
    alf = -math.atan(vC3/np.linalg.norm(vp))

    #Find some unit vectors.
    #vhat is the unit vector in the v direction
    #vphat is the equivalent to the x-body unit vector expressed in
    #inertial coordinates
    #ulat is the cross product of C3 and vphat
    vhat = v/np.linalg.norm(v)
    vphat = vp/np.linalg.norm(vp)
    ulat = np.cross(C3,vphat)

    #Relabel the unit vectors that also happen to be the
    #body unit vectors. These new body unit vectors are labeled
    #xbhat, ybhat, and zbhat respectively and don't necessarily line
    #up with the lab fram unit vectors
    xbhat = vphat
    ybhat = ulat
    zbhat = C3

    #In this brief section the forces are calculated
    #Calculate the lift force
    Flift = CL(alf)*A*rho*np.dot(v,v)/2*np.cross(vhat,ybhat)
    #Drag force
    Fdrag = CD(alf)*A*rho*np.dot(v,v)/2*(-vhat)
    #Gravitational force
    #Note: this is the only force known a priori in the lab frame
    Fgrav = np.array([0,0,-m*g])
    #Total force
    Ftotal = Flift + Fdrag + Fgrav

    #Find the angular velocities in the body frame
    wInC = np.array([phiDot*math.cos(tht),thtDot,phiDot*math.sin(tht)+gamDot])
    #Then in the lab frame
    wInN = np.dot(np.transpose(Tnc),wInC)

    #Find the components of the angular velocities (as expressed in the lab
    #frame) along each of the body framy unit vectors
    wxb = np.dot(wInN,xbhat)
    wyb = np.dot(wInN,ybhat)
    wzb = np.dot(wInN,zbhat)

    #Calculate the spin parameter. This is only used if
    #we can either calculate the Robins-Magnus force or we have
    #more accurate data about drag coefficients
    spinParam = wzb*d/(2*np.linalg.norm(v))

    #In this section we calculate each torque (moment)
    #These torque names are shortened from, e.g, "Tau_x_body" to tauxb
    tauxb = CR(wxb,wzb)*A*rho*d*np.dot(v,v)/2*xbhat
    tauyb = CM(alf,wyb)*A*rho*d*np.dot(v,v)/2*ybhat
    tauzb = CN(wzb)*A*rho*d*np.dot(v,v)/2*np.array([0,0,1])
    #IMPORTANT NOTE: tauxb and tauyb are expressed in the N, since that
    #is the frame of the hat vectors used to construct them. However,
    #The time derivatives of the body angles must be calculated in the
    #body (frisbee) frame. tauzb is constructed so as to already be in
    #the body frame.
    #Now calculate the total torque (moment) in the body frame
    M = np.dot(Tnc,tauxb)+np.dot(Tnc,tauyb)+tauzb

    #Set moments to 0, if you want to
    #M = np.zeros(3)

    #Remove the wind from the velocities
    vx = vx - vwind[0]
    vy = vy - vwind[1]
    vz = vz - vwind[2]

    #Calculate all of the derivatives
    paramsDot[0] = vx # x velocity
    paramsDot[1] = vy # y velocity
    paramsDot[2] = vz # z velocity
    paramsDot[3] = Ftotal[0]/m # x acceleration
    paramsDot[4] = Ftotal[1]/m # y acceleration
    paramsDot[5] = Ftotal[2]/m # z acceleration
    paramsDot[6] = phiDot # phi angular velocity
    paramsDot[7] = thtDot # theta angular velocity
    paramsDot[8] = (M[0]
                    +2*Ixy*thtDot*phiDot*math.sin(tht)
                    -Iz*thtDot*(phiDot*math.sin(tht)+gamDot)
                    )/(Ixy*math.cos(tht)) #phi angular acceleration
    paramsDot[9] = (M[1]
                    +Iz*phiDot*math.cos(tht)*(phiDot*math.sin(tht)+gamDot)
                    -Ixy*phiDot*phiDot*math.cos(tht)*math.sin(tht)
                    )/Ixy # theta angular acceleration
    paramsDot[10] = (M[2]
                     -Iz*paramsDot[8]*math.sin(tht)
                     -Ixy*thtDot*phiDot*math.cos(tht)
                     )/Iz # gamma angular acceleration
    paramsDot[11] = gamDot # gamma angular velocity

    #Return the derivatives
    #if alf>pi/6:
    #    print(alf*180/pi,phi*180/pi,tht*180/pi)
    return paramsDot


####################################################################


#Below are functions used to animate the flight
#Pay no attention to them
def updateTrajectory(num,trajLines,lines):
    for line, data in zip(lines,trajLines):
        line.set_data(data[0:2,:num])
        line.set_3d_properties(data[2,:num])
    return line

####################################################################


def plotAngles(solution,t,nsteps):
    fig = plt.figure()
    plt.subplots_adjust(hspace=0.00001)

    ax1 = plt.subplot(211)
    ax1.plot(t,solution[:,6]*180.0/pi)
    ax1.set_ylabel("Left   Right",fontsize=28)
    ax1.set_xticklabels([])
    ax1.set_ylim(-np.max(np.fabs(solution[:,6]))*180.0/pi,np.max(np.fabs(solution[:,6]))*180.0/pi)

    ax2 = plt.subplot(212)#,sharex=ax1)
    ax2.plot(t,solution[:,7]*180.0/pi)
    ax2.set_ylabel("Pitch Angle",fontsize=28)
    fig.tight_layout()
    fig.savefig(throwname+"angles.png")
    plt.show()

####################################################################


def plotTrajectory(solution,nsteps):
    #Make some shiny 3D parametric plots of x, y, and z
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()

    ax = fig.gca(projection='3d')
    max_dim = max(np.array([max(np.fabs(solution[:,0])),
                            max(np.fabs(solution[:,1]))]))
    ax.set_xlim3d(0,max_dim)
    ax.set_ylim3d(-max_dim/2,max_dim/2)
    ax.set_zlim3d(0,max_dim/2)#max(solution[:,2]))

    fig.tight_layout()
    #Use the following two lines to plot the trajectory immediately
    #full_line = ax.plot(solution[:,0],solution[:,1],
    #                    solution[:,2],label='Flight Path',color='b')
    #ax.legend()#doesn't work unless we draw the static trajectory as well

    #Use the following lines to animate the trajectory
    #lines = [ax.plot(solution[:,0],solution[:,1],solution[:,2])[0]]
    lines = [ax.plot(solution[:,0],solution[:,1],solution[:,2])[0],
             ax.plot(solution[:,0],solution[:,1],solution[:,2]*0,linestyle='--')[0],
             ax.plot(solution[:,0],solution[:,1]*0+max_dim/2,solution[:,2],linestyle=':')[0]]
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    
    #trajectory = np.array([[solution[:,0],solution[:,1],solution[:,2]]])
    trajectory = np.array([[solution[:,0],solution[:,1],solution[:,2]],
                           [solution[:,0],solution[:,1],solution[:,2]*0],
                           [solution[:,0],solution[:,1]*0+max_dim/2,solution[:,2]]])
    lines_ani = animation.FuncAnimation(fig, updateTrajectory,
                                       nsteps,fargs=(trajectory,lines),
                                       interval=0.0001,blit=False)
    print "nsteps = ",nsteps
    #lines_ani.save(throwname+'throw.gif',fps=20)
    plt.show()


####################################################################
    

if __name__ == '__main__':
    main()#This is used so that the program always jumps to main()
