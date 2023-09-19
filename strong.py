import numpy as np
import matplotlib.pyplot  as plt
from sys import argv

time, nodes = np.loadtxt("strong.data", unpack=True)

N = argv[1]
Np = int(argv[2]) + 1
speedup = np.zeros(Np)
for i in range(Np):
    speedup[i]=(time[0]/time[i])

parallel = np.zeros(Np)
for i in range(Np):
    parallel[i]=(time[0]/time[i])/nodes[i]

plt.plot(np.linspace(1,100, 10),np.linspace(1,100, 10),label="Theorical",linestyle="dashed", color="k")
plt.scatter(nodes,speedup, marker='o', color="r", facecolor='none', label="Experimental")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# Processes")
plt.ylabel("speedup")
plt.xlim(1, 100)
plt.legend()
plt.savefig("speedup.png",dpi=200)

plt.cla()

plt.plot([1,100],[1,1],label="Theoretical",linestyle="dashed", color="k")
plt.scatter(nodes,parallel, marker='o', color="r", facecolor='none', label="Experimental")

plt.xlabel("# Processes")
plt.ylabel("Parallel efficiency")
plt.xscale("log")
plt.xlim(1, 100)
plt.legend()
plt.savefig("parallel.png",dpi=200)