Ns = 2000
rep = 10

all: scaling

random.txt: random.cpp
	g++ $< -o random.x;\
	./random.x ${Ns}
scaling: scaling.cpp NBodies.cpp NBodies.h scaling.sh random.txt parallel.py speedup.py
	mpic++ $< NBodies.cpp -o scaling.x;\
	bash scaling.sh ${rep};\
	python3 parallel.py;\
	python3 speedup.py
clean:
	rm -f *.x *.txt *.png
