from numpy import loadtxt, empty, pi, sum, zeros, arange
from numpy import linalg as LA
import matplotlib.pyplot as plt
plt.style.use('dark_background')
from sys import argv
from tqdm import tqdm

N = int(argv[1])    # Number of particles
dt = float(argv[2])    # Size step
jump = int(argv[3])    # Size jump
nP = int(argv[4])

x, y, z, vx, vy, vz = loadtxt("Evolution.txt", unpack=1)    # Evolution

mass = zeros(N)

for ii in range(nP):
    mass[ii*(N//nP):(ii+1)*(N//nP)] = loadtxt(f'data{ii}.txt', usecols = 6, unpack = True)

steps = x.size // N # Printed steps
G = 4*pi**2
t = arange(0, steps*jump*dt, jump*dt)
def TotalEnergy(m, x, y, z, vx, vy, vz):
    global N, G
    U = 0.
    rel = zeros(2)
    v2 = vx**2+vy**2+vz**2
    Ekin = 0.5*sum(mass*v2)
    for ii in range(N):
        for jj in range(ii+1,N):
            rel = [x[ii]-x[jj], y[ii]-y[jj], z[ii]-z[jj]]
            U += -G*(m[jj]*m[ii])/LA.norm(rel)
    return Ekin+U

Energy = empty(steps)
Bar = tqdm(total = steps) #Bar changing
for ii in range(steps):
    Energy[ii] = TotalEnergy(mass, x[ii*N:(ii+1)*N], y[ii*N:(ii+1)*N], z[ii*N:(ii+1)*N], vx[ii*N:(ii+1)*N], vy[ii*N:(ii+1)*N], vz[ii*N:(ii+1)*N])
    Bar.update(1)
Bar.close()
fig, ax = plt.subplots( figsize=(10,7))

ax.plot(t[:len(Energy)], Energy, color='mediumslateblue')
ax.set_title('Energy')
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$Energy$')

plt.show()