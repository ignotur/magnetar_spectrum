
CC = g++
TARGET = magpies.out
LLIST = photon.o monte_carlo.o mathematics.o planck.o main.o
TLIST = photon.o monte_carlo.o mathematics.o planck.o test.o

all: $(LLIST)
		$(CC) -O3 $? -o $(TARGET)
%.o: %.cpp
		g++ -c $< -o $@ 
clean:
		rm *.exe *.o *.out
test: $(TLIST)
		$(CC) $? -o tests.out
		./tests.out
