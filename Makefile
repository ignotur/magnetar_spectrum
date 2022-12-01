all: photon.cpp monte_carlo.cpp mathematics.cpp planck.cpp 
		g++ -O3 photon.cpp monte_carlo.cpp mathematics.cpp planck.cpp -o magpies.out
clean:
		rm *.exe *.o *.out
