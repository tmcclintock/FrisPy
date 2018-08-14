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
#math, numpy, scipy, matplotlib, and emcee
import math
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import emcee
import string

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
######### BELOW HERE IS WHERE YOU CONTROL THE KIND OF THROW ########
####################################################################


#Declare initial conditions
#Note: it is assumed that the initial velocity is in the x-z lab frame plane
#Note: a negative gamDot means the disc is spinning clockwise looking down
#i.e. a backhand release. Thus a positive gamDot is a forehand.
x0,y0,z0 = 0,0,1 #m; (x,y,z) positions. 
speed0 = 12.5    #m/s; initial speed
phi0  = 0.0      #degrees; roll angle (+ means OI BH/IO FH, - is opposite)
tht0 = 0.0       #degrees; polar angle tilt (- means the throw is airbounce-y)
fpAng0 = 0       #degrees; flight path angle between v and x-lab
phiDot0 = 0.0    #rad/s; roll anglular momentum (+ means tilting right)
thtDot0 = 0.0    #rad/s; pitch angular momentum (+ means tilting backwards)
gamDot0 = -50     #rad/s; spin on the disc (- is BH, + is FH)

####################################################################


#Declare a vector that holds the wind
vwind = np.array([0,0,0]) #m/s; (x,y,z) three components of the wind velocity


####################################################################


#Define constant coefficients
#All of these numbers are from Hummell (2003)
CL0 = 0.3331   #Lift coefficient constant 0
CL1 = 1.9124   #Lift coefficient constant 1
CD0 = 0.1769   #Drag coefficient constant 0
CD1 = 0.685    #Drag coefficient constant 1
CR1x = 0.0125  #x-axis Roll angle torque coefficient constant 1
CR1z = 0.00171 #z-axis Roll angle torque coefficient constant 1
CM0 = 0.0821   #Pitch angle torque (moment) coefficient constant 0
CM1a = 0.4338  #Pitch angle y-body torque coefficient constant 1
CM1y = 0.0144  #y-body Pitch angle torque coefficient constant 1
CN1z = 0.0000341 #z-body Spin down torque coefficient constant 1
Cvec0 = np.array([CL0,CL1,CD0,CD1,CR1x,CR1z,CM0,CM1a,CM1y,CN1z])
Cvec = np.copy(Cvec0)

names = [r"C$_{L0}$",r"C$_{L1}$",r"C$_{D0}$",r"C$_{D1}$",r"C$_{\omega_{X}}$",r"C$_{\omega_{Z},X}$",r"C$_{Y0}$",r"C$_{Y1}$",r"C$_{\omega_{Y}}$",r"C$_{\omega_{Z},Z}$"]

####################################################################


#Read in the throw data
#Note
f = open("throw_data.txt",'r')
nlines=0
data = np.zeros((13))
for line in f:
    dataline = [float(s) for s in line.split()]
    data = np.vstack((data,dataline))
    nlines+=1
data = np.delete(data,0,0)
print "Input data has "+str(len(data))+" timesteps of "+str(len(data[0]))+" variables."

#Define the uncertainties in the final observed state
#Note: these numbers are completely made up right now
dxf, dyf, dzf = 0.2, 0.2, 0.01
dvxf, dvyf, dvzf = 0.1, 0.1, 0.1
dphif, dthtf = 0.01, 0.1
dphiDotf, dthtDotf = 0.02, 0.01
dgamDotf, dgamf = 1.0, 5.0
dtimef = 0.1

delta = np.array([dxf,dyf,dzf,dvxf,dvyf,dvzf,dphif,dthtf,dphiDotf,dthtDotf,dgamDotf,dgamf,dtimef])


####################################################################
def main(Cvec,nlines):
    #lnprob(Cvec,data, delta, nlines)
    
    #Here we will do the montecarlo search for the best set of C coefficient values
    #The terminology here follows that for the "emcee" program.
    #See Foreman-Mackey et al. (2012)

    #Define the number of free parameters (C coefficients)
    ndim = len(Cvec)
    
    #First define the number of walkers
    nwalkers = 26

    #Define our initial positions for the walkers
    p0 = np.zeros((nwalkers,len(Cvec)))
    for i in xrange(nwalkers):
        for j in xrange(len(Cvec)):
            p0[i,j] = Cvec[j]*(0.5+np.fabs(np.random.rand()))

    #Initialize the sampler with the chosen specs
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[data, delta, nlines])

    #Run nburnin trials as a burn-in.
    nburnin = 50
    pos, prob, state = sampler.run_mcmc(p0,nburnin)

    #print("Position after burn in:",pos)

    #Reset the chain to remove the burn-in samples
    sampler.reset()

    #Starting from the final position in the burn-in chain,
    #sample for 1000 trials
    ntrials = 200
    sampler.run_mcmc(pos,ntrials,rstate0=state)

    #Print out the mean acceptance fraction
    print("Mean acceptance fraction:",np.mean(sampler.acceptance_fraction))

    #Estimate the integrated autocorrelation time for the time series in each
    #parameter
    print("Autocorrelation time:",sampler.get_autocorr_time())

    #makehistograms(sampler,Cvec)
    #makecontours(sampler)
    makecontoursPartial(sampler,4)
    
def makecontours(sampler):
    nvars = len(sampler.flatchain[0,:])
    ncols = nvars
    nrows = nvars
    fig,axes = plt.subplots(nrows,ncols)

    #Delete unused spots in the subplot
    for i in range(0,nrows):            
        for j in range(i+1,nrows):
            fig.delaxes(axes[i][j])

    #Create contour and histograms
    for i in range(0,nrows):            
        for j in range(0,i+1):
            xdata = sampler.flatchain[:,i]
            ydata = sampler.flatchain[:,j]
            if i==j:
                axes[i][j].hist(xdata,100,range=(0,Cvec0[i]*3))
            else:
                h2d,xe,ye = np.histogram2d(xdata,ydata)
                x = (xe[:-1]+xe[1:])/2.
                y = (ye[:-1]+ye[1:])/2.
                axes[i][j].contour(x,y,h2d.T,3)
                axes[i][j].set_xlim(left=0,right=3.0*Cvec0[j])
                axes[i][j].set_ylim(bottom=0,top=3.0*Cvec0[i])
            if i==nrows-1:
                axes[i][j].set_xlabel(names[j])
            if j==0:
                axes[i][j].set_ylabel(names[i])
            axes[i][j].ticklabel_format(style='sci',axis='x')
            axes[i][j].xaxis.get_major_formatter().set_powerlimits((0,1))
    fig.tight_layout()
    plt.savefig("contourFull.png")
    #plt.show()
    
def makecontoursPartial(sampler,length):
    nvars = length
    ncols = nvars
    nrows = nvars
    fig,axes = plt.subplots(nrows,ncols)

    #Delete unused spots in the subplot
    for i in range(0,nrows):            
        for j in range(i+1,nrows):
            fig.delaxes(axes[i][j])

    #Create contour and histograms
    for i in range(0,nrows):            
        for j in range(0,i+1):
            xdata = sampler.flatchain[:,i]
            ydata = sampler.flatchain[:,j]
            if i==j:
                axes[i][j].hist(xdata,100,range=(0,Cvec0[i]*3))
            else:
                h2d,xe,ye = np.histogram2d(xdata,ydata)
                x = (xe[:-1]+xe[1:])/2.
                y = (ye[:-1]+ye[1:])/2.
                axes[i][j].contour(x,y,h2d.T,3)
                axes[i][j].set_xlim(left=0,right=3.0*Cvec0[j])
                axes[i][j].set_ylim(bottom=0,top=3.0*Cvec0[i])
            if i==nrows-1:
                axes[i][j].set_xlabel(names[j])
            if j==0:
                axes[i][j].set_ylabel(names[i])
            axes[i][j].ticklabel_format(style='sci',axis='x')
            axes[i][j].xaxis.get_major_formatter().set_powerlimits((0,1))
    fig.tight_layout()
    plt.savefig("contourPartial.png")
    plt.show()

def makehistograms(sampler,Cvec):

    FS = 20#fontsize
    Y = 1.04
    fig, axes = plt.subplots(3,5)
    fig.delaxes(axes[2][0])
    fig.delaxes(axes[2][1])
    fig.delaxes(axes[2][2])
    fig.delaxes(axes[2][4])
    fig.delaxes(axes[1][4])
    #CLift
    axes[0][0].hist(sampler.flatchain[:,0],bins=100,range=(0,Cvec0[0]*3))
    axes[0][0].set_title(names[0],fontsize=FS,y=Y)
    axes[1][0].hist(sampler.flatchain[:,1],bins=100,range=(0,Cvec0[1]*3))
    axes[1][0].set_title(names[1],fontsize=FS,y=Y)

    #CDrag
    axes[0][1].hist(sampler.flatchain[:,2],bins=100,range=(0,Cvec0[2]*3))
    axes[0][1].set_title(names[2],fontsize=FS,y=Y)
    axes[1][1].hist(sampler.flatchain[:,3],bins=100,range=(0,Cvec0[3]*3))
    axes[1][1].set_title(names[3],fontsize=FS,y=Y)

    #CX
    axes[0][2].hist(sampler.flatchain[:,4],bins=100,range=(0,Cvec0[4]*3))
    axes[0][2].set_title(names[4],fontsize=FS,y=Y)
    axes[0][2].ticklabel_format(style='sci',axis='x')
    axes[0][2].xaxis.get_major_formatter().set_powerlimits((0,1))
    axes[1][2].hist(sampler.flatchain[:,5],bins=100,range=(0,Cvec0[5]*3))
    axes[1][2].set_title(names[5],fontsize=FS,y=Y)
    axes[1][2].xaxis.get_major_formatter().set_powerlimits((0,1))
    axes[1][2].ticklabel_format(style='sci',axis='x')
    
    #CY
    axes[0][3].hist(sampler.flatchain[:,6],bins=100,range=(0,Cvec0[6]*3))
    axes[0][3].set_title(names[6],fontsize=FS,y=Y)
    axes[1][3].hist(sampler.flatchain[:,7],bins=100,range=(0,Cvec0[7]*3))
    axes[1][3].set_title(names[7],fontsize=FS,y=Y)
    axes[2][3].hist(sampler.flatchain[:,8],bins=100,range=(0,Cvec0[8]*3))
    axes[2][3].set_title(names[8],fontsize=FS,y=Y)
    axes[2][3].ticklabel_format(style='sci',axis='x')
    axes[2][3].xaxis.get_major_formatter().set_powerlimits((0,1))

    #CZ
    axes[0][4].hist(sampler.flatchain[:,9],bins=100,range=(0,Cvec0[9]*3))
    axes[0][4].set_title(names[9],fontsize=FS,y=Y)
    axes[0][4].xaxis.get_major_formatter().set_powerlimits((0,1))
    axes[0][4].ticklabel_format(style='sci',axis='x')
    
    #for i in xrange(len(Cvec)):
    #    plt.hist(sampler.flatchain[:,i],bins=100,range=(0,Cvec0[i]*3))
    #    plt.title("test")
    #    plt.savefig("Histograms/"+names[i]+".png")
    #    #plt.show()
    #    plt.clf()
    #fig.tight_layout()
    plt.subplots_adjust(wspace=0.3,hspace=.3)
    plt.savefig("Histograms.png")
    plt.show()
    print("Generating histograms complete.")

####################################################################


def lnprob(Cvec, data, delta, nlines):
    #Make sure our coefficients make sense (haven't changed signs)
    for i in range(0,len(Cvec)):
        if Cvec[i] < 0:
            return -1e8 #this is a very bad guess
    
    #Calculate the solution at all times of the equations of motion
    solution, times = simulate(Cvec,nlines)

    #Create a vector of the final state solution
    vec_final_calc = solution[nlines-1,:] #vector of calculated final positions and velocities
    t_final = times[nlines-1]
    vec_final_calc = np.append(vec_final_calc,t_final)

    #Calculate the maximum likelihood (-1/2)*chisquared
    #print "solution lengths: ",len(solution),len(solution[0])
    #print "data lengths: ",len(data),len(data[0])
    like = -(1./2.)*np.sum(((solution-data)**2/np.fabs(delta)))

    if math.isnan(like):
        return -1e8

    #Return the likelihood
    return like

####################################################################

#Define our main funciton
def simulate(Cvec,nlines):

    #Set the initial conditions to be fed to the integrator
    x = x0  #m; x-lab displacement
    y = y0  #m; y-lab displacement
    z = z0  #m; height above the ground
    vx = speed0*math.cos(fpAng0*2*pi/360) #initial x-component velocity
    vy = 0.                               #initial y-component velocity
    vz = speed0*math.sin(fpAng0*2*pi/360) #initial z-component velocity
    phi = phi0*2*pi/360     #degrees; roll angle 
    tht = tht0*2*pi/360     #degrees; pitch angle 
    phiDot = phiDot0        #radians/s; roll angular momentum
    thtDot = thtDot0        #radians/s; pitch anglular momentum
    gamDot = gamDot0        #radians/s; initial z-body angular momentum (spin)
    gam = 0.                #radians; spin angle (has no effect on computation)
    
    #Construct the initial parameter array
    params = np.array([x,y,z,vx,vy,vz,phi,tht,phiDot,thtDot,gamDot,gam])

    #Declare the initial and final times
    ti = 0.0   #s; initial time
    tf = 3.5  #s; final time

    #Declare the number of time steps and find the step size
    nsteps = nlines
    dt = (tf-ti)/nsteps

    #Declare the array that contains the times
    #print("Timestep = %4.3f seconds" % dt)
    t = np.arange(ti,tf,dt)
    times = [[time] for time in t]

    #Integrate the solution
    #solution is an nsteps X 12 size array containing the solution at all steps
    solution = odeint(eqOfMotion,params,t,args=(Cvec,))

    #Append the time values at the end of the solution array
    solution = np.append(solution,times,1)

    #Return the solution and time
    return solution, t

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
def CLa(a,CL1):
    #Lift force coefficient
    #function of angle of attack
    return CL1*a#math.sin(a) #from Hummell (2003)

def CL(a,CL0,CL1):
    #Total lift force coefficient.
    return CL0 + CLa(a,CL1)


def CDa(a,CD1):
    #Drag force coefficient
    #function of angle of attack
    a0 = pi/45
    #According to Hummell (2003) this is the angle for which minimal drag is
    #achieved. This actually comes from data.
    return CD1*(a-a0)*(a-a0) #Hummell (2003)

def CD(a,CD0,CD1):
    #Total drag for coefficient
    return CD0 + CDa(a,CD1)

def CRx(wx,CR1x):
    #Torque in the x-body direction (roll moment)
    #functoon of x-body angular velocity (roll angle)
    return -CR1x*wx #Hummell (2003) long flight

def CRz(wz,CR1z):
    #Torque in the x-body direction (roll moment)
    #function of z-body angular velocity (spin)
    return -CR1z*wz  #Hummell (2003)

def CR(wx,wz,CR1x,CR1z):
    #Total x-body torque
    return CRx(wx,CR1x) + CRz(wz,CR1z)

def CMa(a,CM1a):
    #Torque in the y-body direction (pitch moment)
    #function of the angle of attack
    return CM1a*math.sin(a)#really not a good parameterization

def CMy(wy,CM1y):
    #Torque in the y-body direction (pitch moment)
    #function of the pitch angle angular velocity
    return -CM1y*wy #Hummell (2003) long flight

def CM(a,wy,CM0,CM1a,CM1y):
    #Total y-body torque
    return CM0 + CMa(a,CM1a) + CMy(wy,CM1y)

def CNz(wz,CN1z):
    #Torque in the z-body direction (spin down moment)
    #function of the yaw angle angular velocity (spin about the z-body axis)
    return -CN1z*wz #Hummell (2003) long flight

def CN(wz,CN1z):
    #Total z-body torque. Just an alias for CNz made for consistent notation
    return CNz(wz,CN1z)

####################################################################


#Define the differential equation for each variable
def eqOfMotion(params,t,Cvec):
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
    #TOM: I THINK that I want -atan, but it might be +atan
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
    Flift = CL(alf,CL0,CL1)*A*rho*np.dot(v,v)/2*np.cross(vhat,ybhat)
    #Drag force
    Fdrag = CD(alf,CD0,CD1)*A*rho*np.dot(v,v)/2*(-vhat)
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
    tauxb = CR(wxb,wzb,CR1x,CR1z)*A*rho*d*np.dot(v,v)/2*xbhat
    tauyb = CM(alf,wyb,CM0,CM1a,CM1y)*A*rho*d*np.dot(v,v)/2*ybhat
    tauzb = CN(wzb,CN1z)*A*rho*d*np.dot(v,v)/2*np.array([0,0,1])
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
    return paramsDot


####################################################################


#Below are functions used to animate the flight
#Pay no attention to them
def update_trajectory(num,traj,line):
    line.set_data(traj[0:2,:num])
    line.set_3d_properties(traj[2,:num])
    return line

####################################################################


def makePlots(solution,nsteps):
    #Make some shiny 3D parametric plots of x, y, and z
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    max_dim = max(np.array([max(np.fabs(solution[:,0])),
                            max(np.fabs(solution[:,1]))]))
    ax.set_xlim3d(0,max_dim)
    ax.set_ylim3d(-max_dim/2,max_dim/2)
    ax.set_zlim3d(0,max_dim/2)#max(solution[:,2]))

    #Use the following two lines to plot the trajectory immediately
    full_line = ax.plot(solution[:,0],solution[:,1],
                        solution[:,2],label='Flight Path',color='b')
    ax.legend()#doesn't work unless we draw the static trajectory as well

    #Use the following lines to animate the trajectory
    line = ax.plot(solution[:,0],solution[:,1],solution[:,2])[0]
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    

    trajectory = np.array([solution[:,0],solution[:,1],solution[:,2]])
    #line_ani = animation.FuncAnimation(fig, update_trajectory,
    #                                   nsteps,fargs=(trajectory,line),
    #                                   interval=0.0001,blit=False)

    plt.show()


####################################################################
    

if __name__ == '__main__':
    main(Cvec,nlines)#This is used so that the program always jumps to main()
