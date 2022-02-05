import numpy as np
import matplotlib.pyplot  as plt
import sys

Np = sys.argv[1]

x = np.loadtxt("strong1.txt")[:,2]
y1 = np.loadtxt("strong1.txt")[:,1]
y2 = np.loadtxt("strong.txt")[:,1]

plt.scatter(x,y1/x,marker='.',color="r",label="1 Process ")
plt.scatter(x,y2/x,marker='.',color="k",label= Np + " Processes ")

plt.xlabel("Number of particles")
plt.ylabel("Time [s]")
plt.title("Strong Scaling")
plt.legend()
plt.savefig("strong.png",dpi=200)
