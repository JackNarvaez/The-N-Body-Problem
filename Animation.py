from turtle import color
from numpy import loadtxt, empty, max
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from sys import argv

N = int(argv[1])    # Number of particles
dt = float(argv[2])    # Size step
jump = int(argv[3])    # Size jump

x, y, z = loadtxt("Evolution.txt", unpack=1)    # Evolution
steps = x.size // N # Printed steps

System = empty((N, 3, steps))

for jj in range(0,N):
    for ii in range(0,steps):
        System[jj,0, ii]= x[N*ii+jj]
        System[jj,1, ii]= y[N*ii+jj]
        System[jj,2, ii]= z[N*ii+jj]

boundaries = max(System)

def Gen_orbit(particle, time, dims=3):
    """-----------------------------------------------
    Create the orbit line using Evolution.txt
    --------------------------------------------------
    Arguments:
    particles   :  Id of particle
    time        :  Number of time steps.
    dims        :  Number of dimensions the orbit has.
    -----------------------------------------------"""
    lineData = empty((dims, time))
    lineData[:, 0] = System[particle, :,0]
    for index in range(1, steps):
        lineData[:, index] = System[particle, :, index]

    return lineData


def update_orbit(num, dataOrbits, orbits):
    for orbit, data in zip(orbits, dataOrbits):
        message.set_text(f"Time = {num*jump*dt: .2f} years")
        orbit.set_data(data[0:2, :num])
        orbit.set_3d_properties(data[2, :num])
    return orbits

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

message = ax.text2D(0.00, 0.9, "", transform=ax.transAxes)
# Nth orbits
data = [Gen_orbit(ii, steps, 3) for ii in range(N)]

# Creating N orbit objects.
orbits = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
ax.plot(0, 0, "o", color="k")

# Setting the axes properties
ax.set_xlim3d([-boundaries, boundaries])
ax.set_xlabel('x [au]')
ax.set_ylim3d([-boundaries, boundaries])
ax.set_ylabel('y [au]')
ax.set_zlim3d([-boundaries, boundaries])
ax.set_zlabel('z [au]')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_orbit, steps, fargs=(data, orbits),
                                   interval=1, blit=False)

line_ani.save('Evolution.gif', writer='pillow', fps=30)
plt.show()