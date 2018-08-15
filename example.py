"""
This file contains an example of how to use FrisPy to visualize a simulated trajectory.
"""
import FrisPy

#Some example initial conditions are in simple_initial_conditions.txt, located a the top level of this repository.
disc = FrisPy.create_disc(filename = "simple_initial_conditions.txt")
times, trajectory = FrisPy.get_trajectory(disc)

#Try plotting
try:
    import matplotlib.pyplot as plt
    x,y,z = trajectory.T
    plt.plot(x,z)
    plt.xlabel("x [meters]")
    plt.ylabel("z [meters]")
    plt.show()
except Exception:
    print "Matplotplot not installed. Cannot visualize example."
