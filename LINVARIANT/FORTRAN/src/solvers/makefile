include makefile.include
 
exe = LINVARIANT.x
 
CORE_DIR   = ../core
IO_DIR     = ../io

OBJ=$(%.f90=.o)

.PHONY : all linvariant io

all : mc.o md.o pt.o

%.o : %.f90 linvariant io
	$(FC) -c $(FFLAGS) $< -o $@

linvariant:
	rsync -ru $(CORE_DIR)/*.mod .

io:
	rsync -ru $(IO_DIR)/*.mod .
