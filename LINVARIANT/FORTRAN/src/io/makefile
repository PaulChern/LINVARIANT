include makefile.include
 
exe = LINVARIANT.x
 
CORE_DIR   = ../core

OBJ=$(%.f90=.o)

.PHONY : all linvariant

all : Inputs.o Outputs.o

Outputs.o : Outputs.f90 linvariant
	$(FC) -c $(FFLAGS) $< -o $@

Inputs.o : Inputs.f90 FileParser.o linvariant
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

linvariant:
	rsync -ru $(CORE_DIR)/*.mod .
