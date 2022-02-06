import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import mpl_toolkits.mplot3d as Axes3D
import sys

N = int(sys.argv[1])

x = np.loadtxt("Evolution.txt")[:,0]
y = np.loadtxt("Evolution.txt")[:,1]
z = np.loadtxt("Evolution.txt")[:,2]
# Limits for the plot
#print(x, y, z)
fig = plt.figure(figsize=(10, 7))
ax = fig.gca(projection='3d')
ax.scatter3D(x, y, z,color='white')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#ax.set_xlim3d(-boundary, boundary)
#ax.set_ylim3d(-boundary, boundary)
#ax.set_zlim3d(-boundary, boundary)
#plt.savefig(path+'NBody-output.jpg')
plt.show()
