#Set Compiler and Flags
#---------------------------------------------------------------
FC90 = mpif90
FLAGS = -O2
LIBS = -l lapack
#Set directory
#---------------------------------------------------------------
SRCDIR = ./src
BINDIR = ./bin
OBJDIR = ./obj
MODDIR = ./mod
#Set dependencies
#---------------------------------------------------------------
MOD = $(wildcard $(SRCDIR)/mod_*.F90)
SRC = $(wildcard $(SRCDIR)/src_*.F90)

CONS = $(patsubst $(SRCDIR)/%.F90,./%.mod,$(CON))
BINS = $(BINDIR)/PIC_exe
MODS = $(patsubst $(SRCDIR)/%.F90,./%.mod,$(MOD))

all: $(BINS)

$(BINS): $(MOD) $(SRC)
	$(FC90) $(FLAGS) $(MOD) $(SRC) -o $@ $(LIBS)

#%.mod: $(CON)
#	$(FC90) $(FLAGS) $(CON) -o $@

#$(MODDIR)/%.mod: $(SRCDIR)/%.F90
#	$(FC90) $(FLAGS) -c $< -o $@
#	rm *.mod

clean:
	rm $(BINS) $(MODS)

play:
	echo $(MODS)
