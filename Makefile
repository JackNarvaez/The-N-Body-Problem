Ns = 1000 		#Number of particles in scailing
Ng = 300		#Number of particles in the galaxy
i = 0.1 		#Inclination of the galaxy (rads/pi)
w = 0.6 		#Angle in the xy plane of the galaxy (rads/pi)
rep = 100 		#Repetitions in the scaling
Npv =$(shell nproc) 	#Maximun number of threads
Np = $(shell echo $$(($(Npv) / 2 ))) #Maximun number of cores
steps = 70000		#Steps in galaxy animation
dt = 0.0001		#Step size in galaxy animation
jump = 500		#Every jump steps it saves a frame
rad = 5000		#Radius of Galaxy (AU)

all: Galaxy
Random: random.cpp
	g++ $< -o random.x;\
	./random.x ${Ns}

Galaxy: Evolution.cpp NBodies.cpp NBodies.h galaxy.py Animation.py
	python3 galaxy.py ${Ng} ${i} ${w} ${rad} ${Np};\
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} ${Ng};\
	python3 Animation.py ${Ng} ${dt} ${jump} ${rad} 

SagA: Evolution.cpp NBodies.cpp NBodies.h SagA.data Animation.py
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x 500000 0.0001 500 SagA.data;\
	python3 Animation.py 14 0.0001 500 0 	#Default parameters for SagA

Scaling: scaling.cpp NBodies.cpp NBodies.h scaling.sh Random parallel.py speedup.py
	mpic++ $< NBodies.cpp -o scaling.x;\
	bash scaling.sh ${rep};\
	python3 parallel.py ${Ns};\
	python3 speedup.py ${Ns}

Strong: strong_scaling.sh scaling.cpp random.cpp NBodies.cpp NBodies.h fit.cpp
	mpic++ scaling.cpp  NBodies.cpp -o scaling.x
	bash $<  ${Np};\
	bash $< 1;\
	g++ -std=c++17  fit.cpp -lgsl -lgslcblas -o fit.x;\
	echo Number of processes: ${Np};\
	./fit.x ${Np};\
	echo Number of processes: 1;\
	./fit.x 1;\
	python3 strong.py ${Np}

clean:
	rm -f *.x *.txt *.png *.out *.gif
