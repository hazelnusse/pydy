all: bicycle.o

bicycle.o : bicycle.c bicycle.h
	        gcc -Wall -O3 -funroll-loops -c bicycle.c

clean:
	        rm -f bicycle.o
