N = 20
rep = 10
Np = 4
steps = 200000
dt = 0.0002
jump = 100

all: scaling

Random.txt: random.cpp
	g++ $< -o random.x;\
	./random.x ${N}

Galaxy: Evolution.cpp NBodies.cpp NBodies.h galaxy.py
	python3 galaxy.py 4 0 0 > Galaxy.txt;\
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} Galaxy.data

SagA: Evolution.cpp NBodies.cpp NBodies.h SagA.data
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} SagA.data

EvolveRandom: Evolution.cpp NBodies.cpp NBodies.h Random.txt
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} Random.txt 

scaling: scaling.cpp NBodies.cpp NBodies.h scaling.sh Random.txt parallel.py speedup.py
	mpic++ $< NBodies.cpp -o scaling.x;\
	bash scaling.sh ${rep};\
	python3 parallel.py;\
	python3 speedup.py

strong: strong_scaling.sh scaling.sh random.cpp
	mpic++ scaling.cpp  NBodies.cpp -o scaling.x
	bash $<  ${Np};\
	bash $< 1;\
	python3 strong.py Np

clean:
	rm -f *.x *.txt *.png *.out
