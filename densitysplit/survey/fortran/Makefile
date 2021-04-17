# Compiler choice
FC = gfortran
IFC = ifort
MPIFC = gfortran

MLIBS = -lm

#################

EXES = gaussian_filter.exe tophat_filter.exe

all: $(EXES)

gaussian_filter.exe: gaussian_filter.o
	$(FC) gaussian_filter.o -o gaussian_filter.exe $(MLIBS)

tophat_filter.exe: tophat_filter.f90
	$(FC) -fopenmp tophat_filter.f90 -o tophat_filter.exe $(MLIBS)


gaussian_filter.o: gaussian_filter.f90
	$(FC) -c gaussian_filter.f90


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
