LINKFLAGS_FOR = -pedantic -O1 -march=native -fbounds-check
COMP_FOR = gfortran
install:
	${COMP_FOR} ${LINKFLAGS_FOR} farthrest_nodes_subtitutions.f90 -c
	${COMP_FOR} ${LINKFLAGS_FOR} ising_frameworks.f90 -c
	${COMP_FOR} ${LINKFLAGS_FOR} ising_frameworks.o farthrest_nodes_subtitutions.o -o ising_frameworks
all:
	make install
	make execute
	make clean
execute:
	time ./ising_frameworks
clean:;         @rm -f *.o
