Ns = 2000
rep = 10

all: scaling.jpg

File.txt: File.cpp
	g++ $< -o random.x;\
	./random.x ${Ns}
scaling: main.cpp NBodies.cpp NBodies.h scaling.sh File.txt parallel.py speedup.py
	mpic++ $< NBodies.cpp -o scaling.x;\
	bash scaling.sh ${Ns} ${rep};\
	python3 parallel.py;\
	python3 speedup.py
