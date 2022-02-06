Ns = 100
Ng = 10
rep = 10
Npv =$(shell nproc)
Np = $(shell echo $$(($(Npv) / 2 )))
steps = 200000
dt = 0.0002
jump = 100

all: scaling

Random: random.cpp
	g++ $< -o random.x;\
	./random.x ${Ns}

Galaxy: Evolution.cpp NBodies.cpp NBodies.h galaxy.py
	python3 galaxy.py ${Ng} 0 0 > Galaxy.txt;\
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} Galaxy.txt

SagA: Evolution.cpp NBodies.cpp NBodies.h SagA.data
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} SagA.data

EvolveRandom: Evolution.cpp NBodies.cpp NBodies.h Random
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x ${steps} ${dt} ${jump} Random.txt 

scaling: scaling.cpp NBodies.cpp NBodies.h scaling.sh Random parallel.py speedup.py
	mpic++ $< NBodies.cpp -o scaling.x;\
	bash scaling.sh ${rep};\
	python3 parallel.py ${Ns};\
	python3 speedup.py ${Ns}

strong: strong_scaling.sh scaling.sh random.cpp
	mpic++ scaling.cpp  NBodies.cpp -o scaling.x
	bash $<  ${Np};\
	bash $< 1;\
	python3 strong.py ${Np}

clean:
	rm -f *.x *.txt *.png *.out
