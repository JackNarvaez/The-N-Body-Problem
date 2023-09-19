import numpy as np
import matplotlib.pyplot  as plt
from sys import argv

time1, nodes1, time2, nodes2, time3, nodes3 = np.loadtxt("weak.data", unpack=True)

Np = argv[1]
N = len(time1)
strong1 = np.zeros(N)
strong2 = np.zeros(N)
strong3 = np.zeros(N)
bodies = 10*np.arange(1, int(argv[2])+1)
for i in range(N):
    strong1[i]=(time1[i]/time1[0])
    strong2[i]=(time2[i]/time2[0])
    strong3[i]=(time3[i]/time3[0])

plt.scatter(bodies,strong1, marker='o', color="r", facecolor='none', label="10 Proc.")
plt.scatter(bodies,strong2, marker='o', color="b", facecolor='none', label="20 Proc.")
plt.scatter(bodies,strong3, marker='o', color="g", facecolor='none', label="40 Proc.")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# Bodies")
plt.ylabel("Time")
plt.legend()
plt.savefig("weak.png",dpi=200)