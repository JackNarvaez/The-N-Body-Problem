Ns = 2000
rep = 10
Np = 2

all: scaling

Random.txt: random.cpp
	g++ $< -o random.x;\
	./random.x ${Ns}

scaling: scaling.cpp NBodies.cpp NBodies.h scaling.sh Random.txt parallel.py speedup.py
	mpic++ $< NBodies.cpp -o scaling.x;\
	bash scaling.sh ${rep};\
	python3 parallel.py;\
	python3 speedup.py

strong: strong_scaling.sh scaling.sh random.cpp
	mpic++ scaling.cpp  NBodies.cpp -o scaling.x
	bash $<;\
	python3 strong.py

evolution: Evolution.cpp NBodies.cpp NBodies.h Random.txt
	mpic++ $< NBodies.cpp -o Evolution.x;\
	mpirun -np ${Np} ./Evolution.x

clean:
	rm -f *.x *.txt *.png
