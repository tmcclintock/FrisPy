import matplotlib.pyplot as plt

from frispy import Disc

# A negative "theta" is an airbounce, or tilts the disk "back"
# toward the thrower
disc = Disc(vx=10, theta=-0.2)
result = disc.compute_trajectory()
times = result.times
x, y, z = result.x, result.y, result.z

plt.plot(x, z)
plt.show()
