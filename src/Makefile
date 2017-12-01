LINKFLAGS_FOR = -pedantic -O3 -march=native
GNUPLOT_SCRIPT = gp.NMR
COMP_FOR = gfortran
install:
	${COMP_FOR} ${LINKFLAGS_FOR} NMR.f90 -o NMR
all:
	make install
	make execute
	make clean
	make plot
execute:
	time ./NMR -f -g -e -t
plot:
	gnuplot ${GNUPLOT_SCRIPT} 
clean:;         @rm -f *.o *.mod NMR
