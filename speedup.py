from numpy import loadtxt, arange
import matplotlib.pyplot  as plt
import sys

y = loadtxt("scaling.txt")[:,1]
x = arange(1,y.size+1)

N = sys.argv[1]  
T1=y[0]

for i in x:
    y[i-1]=(T1/y[i-1])


plt.plot(x,x,label="Theorical",linestyle="dashed")
plt.scatter(x,y,marker='+',color="k",label="Experimental")

plt.xlabel("# Processes")
plt.ylabel("speedup")
plt.title(f"Speedup with N={N}({y.size//2} cores/{y.size} threads)")
plt.legend()
plt.savefig("speedup.png",dpi=200)