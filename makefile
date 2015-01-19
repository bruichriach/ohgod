FFLAGS=-I/usr/local/include -L/usr/local/lib -cpp -O2 -march=native -mtune=native -ffixed-line-length-none
#FFLAGS= -O3 -m 3 
#FFLAGS=-I/usr/local/include -L/usr/local/lib -cpp -O0 -g -Wall -ffixed-line-length-none -finit-real=nan -finit-integer=666 -fbounds-check

F90=mpif90

default: all

OBJS = params.o global.o optimise.o grid.o system.o operate.o solver.o writeout.o stochastic.o timeav.o


%.o: src/%.F90
	$(F90) $(FFLAGS) -c $< -o src/$@

all: clean folder main

folder:
	mkdir -p 'in'
	mkdir -p 'out'
	mkdir -p 'data'
	mkdir -p 'global'
	mkdir -p 'ierr'


clean:
	rm -f *.mod $(addsuffix .mod, $(basename $(OBJS)))
	rm -f src/*.o $(OBJS) main.o
	rm -f sw.out

main: $(OBJS) main.o
	$(F90) $(FFLAGS) src/*.o -o sw.out

optimise.o: params.o
operate.o: params.o global.o grid.o
grid.o: params.o optimise.o system.o
