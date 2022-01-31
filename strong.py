import numpy as np
import matplotlib.pyplot  as plt

x = np.loadtxt("strong.txt")[:,2]
y = np.loadtxt("strong.txt")[:,1]
    
plt.scatter(x,y,marker='+',color="k",label=f"Experimental")

plt.xlabel("Number of particles")
plt.ylabel("Time [ms]")
plt.title("Strong Scaling (8 processes)")
plt.legend()
plt.savefig("strong.png",dpi=200)
