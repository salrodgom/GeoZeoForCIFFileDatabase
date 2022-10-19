# Salvador R. G. Balestra, 2020
# Makefile for zeo_analyser_angles_OSiO_SiOSi_OOO.f90
FFLAGS = -ffree-form -march=native -O0
PROGRAM = zeo_analyser.exe
install:
	gfortran zeo_analyser.f90 -o $(PROGRAM) $(FFLAGS)
all:
	gfortran zeo_analyser.f90 -o $(PROGRAM) $(FFLAGS)
	./$(PROGRAM)
	rm -rf $(PROGRAM) *.o *.mod
clean:
	rm -rf $(PROGRAM) *.o *.mod
