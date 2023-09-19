Nb = 10000		# Number of bodies in the galaxy
Np = 32     	# Number of processes
i = 0.1 		# Inclination of the galaxy (rads/pi)
w = 0.6 		# Angle in the xy plane of the galaxy (rads/pi)
steps = 10000 	# Evolution steps
jump = 100		# Data storage interval
dt = 0.001		# Time step
rad = 5000		# Radius of Galaxy (AU)

all: PlotGal

Galaxy: Evolution.cpp NBodies.cpp NBodies.h galaxy.py Animation.py
	mpirun -n 4 python3 galaxy.py 100 ${i} ${w} ${rad};\
	mpic++ -std=c++17 $< NBodies.cpp -o Evolution.x;\
	mpirun -np 4 ./Evolution.x 100000 ${dt} ${jump} 100;\
	python3 Animation.py 100 ${dt} ${jump} ${rad}

PlotGal: Animation.py
	python3 $< ${Nb} ${dt} ${jump} ${rad}

SagA: SagA.cpp NBodies.cpp NBodies.h SagA.data Orbits.py
	mpic++ $< NBodies.cpp -o SagA.x;\
	mpirun -np 2 ./SagA.x 100000 0.001 100;\
	python3 Orbits.py 14 0.001 100 0

Strong: strong.py
	python3 $< 1000 6

Weak: weak.py
	python3 $< 10 4

clean:
	rm -f *.x *.txt *.png *.out *.gif
