from numpy import loadtxt
import matplotlib.pyplot  as plt
import sys

Np = sys.argv[1]

x = loadtxt("strong1.txt")[:,2]
y1 = loadtxt("strong1.txt")[:,1]
y2 = loadtxt(f"strong{Np}.txt")[:,1]

plt.scatter(x,y1,marker='.',color="r",label="1 Process ")
plt.scatter(x,y2,marker='.',color="k",label= f"{Np} Processes ")

plt.xlabel("Number of particles")
plt.ylabel("Time [s]")
plt.title("Strong Scaling")
plt.legend()
plt.savefig("strong.png",dpi=200)