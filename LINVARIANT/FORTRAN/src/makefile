.SUFFIXES:
.SUFFIXES: .mod .o .f90
.PRECIOUS: %.f90
# .SECONDEXPANSION:

FUFFIX=.F
SUFFIX=.f90

include makefile.include
include .objects
 
BINDIR = ../../bin
SRCDIR = ../../src
OBJS = $(foreach obj, $(SRCS:.f90=.o), $(obj))

.PHONY: all sources

all: sources
	rm -f $(exe)
	$(MAKE) LINVARIANT
	cp $(exe) $(BINDIR)/$(exe)

LINVARIANT: $(OBJS) main.o
#	$(FC) -o $(exe) $(OBJS) main.o
	$(FCL) -o $(exe) $(OBJS) main.o $(FFT3D) $(LINK) $(LLIBS)

sources:
	rsync -ru $(shell find $(SRCDIR) -name "*.f90") . 


%.o : %.f90 #$(MMA)
	$(FC) -o $@ $(FFLAGS) $(OFLAG) $(INCS) -c $<

#liblinvariant.a: $(linvariant)
#       rm -f liblinvariant.a
#       ar vq liblinvariant.a $(linvariant)

