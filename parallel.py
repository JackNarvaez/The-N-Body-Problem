from numpy import loadtxt, arange
import matplotlib.pyplot  as plt
from sys import argv

y = loadtxt("scaling.txt")[:,1]
x = arange(1,y.size+1)

N = argv[1]  
T1=y[0]

for i in x:
    y[i-1]=(T1/y[i-1])/i

plt.plot([x[0],x[-1]],[1,1],label="Theoretical",linestyle="dashed")
plt.scatter(x,y,marker='+',color="k",label="Experimental")

plt.xlabel("# Processes")
plt.ylabel("Parallel efficiency")
plt.title(f"Parallel Efficiency with N={N} ({y.size//2} cores/{y.size} threads)")
plt.legend()
plt.savefig("parallel.png",dpi=200)