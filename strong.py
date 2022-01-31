import numpy as np
import matplotlib.pyplot  as plt

x = np.loadtxt("strong1.txt")[:,2]
y1 = np.loadtxt("strong1.txt")[:,1]
y2 = np.loadtxt("strong8.txt")[:,1]
 
plt.scatter(x,y1,marker='.',color="r",label=f"1 Process ")
plt.scatter(x,y2,marker='.',color="k",label=f"8 Processes ")

plt.xlabel("Number of particles")
plt.ylabel("Time [s]")
plt.title("Strong Scaling")
plt.legend()
plt.savefig("strong.png",dpi=200)
