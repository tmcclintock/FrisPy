"""
FrisPy is a simple, pure Python implementation of a physical model for a flying disc (frisbee). The default configuration is for a 175 g Ultrastar.
"""

from disc import *

__author__ = "Tom McClintock <tmcclintock89@gmail.com>"
__version__ = "1.0.0"

def create_disc(initial_conditions = None, filename = None, debug=False):
    """Create a disc object and return it to the user.

    Note: either the initial_conditions or a filename must be supplied, but not both.

    Args:
        initial_conditions (array like): Array of initial conditions of all kinematic variables. Optional.
        filename (string): Filename of file containing initial conditions. Optional.
        debug (boolean): Turn the debug option on in the Disc object to print diagnostic information. Optional.

    Returns:
        Disc object: disc object without a model specified.

    """
    if not initial_conditions and not filename:
        raise Exception("Must supply either initial_conditions in an array or the filename that contains them.")
    if initial_conditions and filename:
        raise Exception("Cannot supply both initial_conditions array and a filename for them.")
    if initial_conditions:
        x,y,z,vx,vy,vz,phi,theta,gamma,phidot,thetadot,gammadot = initial_conditions
        disc = Disc(x,y,z,vx,vy,vz,phi,theta,gamma,phidot,thetadot,gammadot,debug)
    else:
        x,y,z,vx,vy,vz,phi,theta,gamma,phidot,thetadot,gammadot = np.loadtxt(filename).T
        disc = Disc(x,y,z,vx,vy,vz,phi,theta,gamma,phidot,thetadot,gammadot,debug)
    return disc

def get_trajectory(disc, times=None, full_trajectory=False):
    """Get the trajectory of the disc object.

    Note: supplying the times and using the default times will result in numerical differences at the floating point level. This is not an issue for the vast majority of applications.

    Args:
        disc (Disc object): the disc object that will be simulated.
        times (array like): array of times in seconds where the trajectory will be evaluated. Optional.
        full_trajectory (boolean): If False, only return an array with x,y,z, otherwise return all kinematic variables. Optional, default is False.

    Returns:
        times (array like): times in seconds where the trajectory is evaluated.
        trajectory (array like): 2D array (N_times x N_variables) containing the kinematic variables at each time.

    """
    if times is None:
        print "No trajectory times supplied. Integrating for 3 seconds with dt=1 ms."
        ti, tf = 0, 3.
        dt = 0.001
    else:
        ti, tf = times[0], times[-1]
        dt = times[1] - times[0]
    if not disc.has_model:
        print "Initializing disc with default model."
        disc.initialize_model()

    #Get the disc's trajectory
    return disc.get_trajectory(ti, tf, dt, full_trajectory)

if __name__ == "__main__":
    d = create_disc(filename="simple_initial_conditions.txt")
    times, trajectory = get_trajectory(d)
    x,y,z = trajectory.T    

    import matplotlib.pyplot as plt
    plt.plot(x,z)
    plt.show()
