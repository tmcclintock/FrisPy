from frispy import Disc
import matplotlib.pyplot as plt

disc = Disc()
result = disc.compute_trajectory()
times = result.times
x, y, z = result.x, result.y, result.z

plt.plot(times, z)
plt.show()
