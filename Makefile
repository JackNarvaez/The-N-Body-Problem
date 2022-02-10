Ns = 100 		#Number of particles in scailing
Ng = 300		#Number of particles in the galaxy
i = 0.1 		#Inclination of the galaxy
w = 0.6 		#Angle in the xy plane of the galaxy
rep = 100 		#Repetitions in the scaling
Npv =$(shell nproc) 	#Number of threads
Np = $(shell echo $$(($(Npv) / 2 ))) #Number of cores
steps = 7000		#Steps in galaxy animation
dt = 0.0001		#Step size in galaxy animation
jump = 40		#Every jump steps it saves a frame

all: Galaxy
Random: random.cpp
	g++ $< -o random.x;\
	./random.x ${Ns}

Galaxy: Evolution.cpp NBodies.cpp NBodies.h galaxy.py Animation.py
	python3 galaxy.py ${Ng} ${i} ${w} > Galaxy.txt;\
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} Galaxy.txt;\
	python3 Animation.py ${Ng} ${dt} ${jump}

SagA: Evolution.cpp NBodies.cpp NBodies.h SagA.data Animation.py
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x 500000 0.0001 500 SagA.data;\
	python3 Animation.py 14 0.0001 500

scaling: scaling.cpp NBodies.cpp NBodies.h scaling.sh Random parallel.py speedup.py
	mpic++ $< NBodies.cpp -o scaling.x;\
	bash scaling.sh ${rep};\
	python3 parallel.py ${Ns};\
	python3 speedup.py ${Ns}

strong: strong_scaling.sh scaling.sh random.cpp NBodies.cpp NBodies.h fit.cpp
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
