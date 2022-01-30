import numpy as np
import matplotlib.pyplot  as plt

#Theory
l1=[1,16]
l2=[1,1]
plt.plot(l1,l2,label=f"Theoretical",linestyle="dashed")

x = np.arange(1,17)
y = np.loadtxt("scaling.txt")[:,1]

T1=y[0]

for i in x:
    y[i-1]=(T1/y[i-1])/i

plt.scatter(x,y,marker='+',color="k",label=f"Experimental")

plt.xlabel("# Processes")
plt.ylabel("Parallel efficiency")
plt.title("Parallel Efficiency with N=2000 (8 cores/16 threads)")
plt.legend()
plt.savefig("parallel.png",dpi=200)
