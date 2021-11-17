import matplotlib.pyplot as plt

from frispy import Disc

disc = Disc(vx=10)
result = disc.compute_trajectory()
times = result.times
x, y, z = result.x, result.y, result.z

plt.plot(x, z)
plt.show()
