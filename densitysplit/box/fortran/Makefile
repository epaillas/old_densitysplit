# Compiler choice
FC = gfortran
IFC = ifort
MPIFC = gfortran

MLIBS = -lm

#################

EXES = tophat_filter.exe

all: $(EXES)

tophat_filter.exe: tophat_filter.f90
	$(FC) -fopenmp procedures.f90 tophat_filter.f90 -o tophat_filter.exe $(MLIBS)


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
