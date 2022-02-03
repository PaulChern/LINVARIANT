include makefile.include
 
exe = LINVARIANT.x
 
SRC_DIR   = ../../../src
BUILD_DIR = ../../build
BIN_DIR   = ../../../bin

linvariant=Constants.o Parameters.o LINVARIANT.o

.PHONY :  EwaldMatrix.o liblinvariant.a

EwaldMatrix.o : EwaldMatrix.f90 liblinvariant.a
	$(FC) -c $(FFLAGS) $< -o $@

liblinvariant.a: $(linvariant)
	rm -f liblinvariant.a
	ar vq liblinvariant.a $(linvariant)

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

