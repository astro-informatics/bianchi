# Makefile for Bianchi library
# Jason McEwen    


# ======== OPTIONS ========

USEPGPLOT = no
#USEPGPLOT = yes


# ======== COMPILER ========

FC      = gfortran
#FC      = f95
#FC      = g95

ifneq ($(USEPGPLOT),yes)
  OPTPGPLOT     = -DNO_PGPLOT
endif
OPT = $(OPTPGPLOT) -DMILLIK -m64


# ======== LINKS ========

PROGDIR      = /Users/jdm/Src

HPIXDIR      = $(PROGDIR)/Healpix
HPIXLIB      = $(HPIXDIR)/lib
HPIXLIBNM    = healpix
HPIXINC      = $(HPIXDIR)/include

S2DIR        = $(PROGDIR)/s2
S2LIB        = $(S2DIR)/lib
S2LIBNM      = s2
S2INC        = $(S2DIR)/include
S2SRC        = $(S2DIR)/src/mod
S2PROG       = $(S2DIR)/src/prog
S2BIN        = $(S2DIR)/bin
S2DOC        = $(S2DIR)/doc

BIANCHIDIR   = $(PROGDIR)/bianchi
BIANCHISRC   = $(BIANCHIDIR)/src/mod
BIANCHIPROG  = $(BIANCHIDIR)/src/prog
BIANCHIINC   = $(BIANCHIDIR)/include
BIANCHIBIN   = $(BIANCHIDIR)/bin
BIANCHIDOC   = $(BIANCHIDIR)/doc
BIANCHILIB   = $(BIANCHIDIR)/lib
BIANCHILIBNM = bianchi

CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
CFITSIOLIBNM = cfitsio

PGPLOTLIB    = $(PROGDIR)/pgplot
PGPLOTLIBNM  = pgplot
X11LIB       = /usr/X11R6/lib
X11LIBNM     = X11


# ======== FFFLAGS ========

FFLAGS  = -I$(HPIXINC) -I$(S2INC)  -I$(BIANCHIINC) -I.


# ======== LDFLAGS ========

ifeq ($(USEPGPLOT),yes)
  LDFLAGSPGPLOT = -L$(PGPLOTLIB) -L$(X11LIB) \
                  -l$(PGPLOTLIBNM) -l$(X11LIBNM)
endif

LDFLAGS =  -L$(BIANCHILIB) -l$(BIANCHILIBNM) \
           -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) \
           -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) \
           $(LDFLAGSPGPLOT)


# ======== PPFLAGS ========

ifeq ($(FC),f95)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -x f95-cpp-input $(OPT)
endif


# ======== OBJECT FILES TO MAKE ========

BIANCHIOBJ = $(BIANCHIINC)/bianchi_sky_mod.o       \
             $(BIANCHIINC)/bianchi_plm1table_mod.o \
             $(BIANCHIINC)/bianchi_error_mod.o


# ======== MAKE RULES ========

default: all

all:     lib prog

lib:	 $(BIANCHILIB)/lib$(BIANCHILIBNM).a 

prog:    $(BIANCHIBIN)/bianchi_sim


$(BIANCHIINC)/%.o: $(BIANCHISRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(BIANCHIINC)

$(BIANCHIINC)/%.o: $(BIANCHIPROG)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 


# Library

$(BIANCHILIB)/lib$(BIANCHILIBNM).a: $(BIANCHIOBJ)
	ar -r $(BIANCHILIB)/lib$(BIANCHILIBNM).a $(BIANCHIOBJ)


# Documentation

docs:
	f90doc_fpp $(BIANCHISRC)/*.f90
	f90doc_fpp $(BIANCHIPROG)/*.f90
	ln_multi $(S2DOC)/s2_*
	ln_multi $(S2DOC)/index_*
	mv *.html $(BIANCHIDOC)/.
	addstyle $(BIANCHIDOC)/bianchi_*

cleandocs:
	rm -f $(BIANCHIDOC)/bianchi_*.html
	rm -f $(BIANCHIDOC)/s2_*.html
	rm -f $(BIANCHIDOC)/index_s2.html


# Cleaning up

clean:	tidy
	rm -f $(BIANCHIINC)/*.mod
	rm -f $(BIANCHIINC)/*.o
	rm -f $(BIANCHILIB)/lib$(BIANCHILIBNM).a
	rm -f $(BIANCHIBIN)/*

tidy:	
	rm -f *.mod
	rm -f $(BIANCHISRC)/*~
	rm -f $(BIANCHIPROG)/*~


# Module dependencies

$(BIANCHIINC)/bianchi_error_mod.o:     $(BIANCHISRC)/bianchi_error_mod.f90
$(BIANCHIINC)/bianchi_plm1table_mod.o: $(BIANCHISRC)/bianchi_plm1table_mod.f90 \
                                         $(BIANCHIINC)/bianchi_error_mod.o
$(BIANCHIINC)/bianchi_sky_mod.o:       $(BIANCHISRC)/bianchi_sky_mod.f90 \
                                         $(BIANCHIINC)/bianchi_plm1table_mod.o \
                                         $(BIANCHIINC)/bianchi_error_mod.o


# Program dependencies and compilation

$(BIANCHIINC)/bianchi_sim.o: $(BIANCHIPROG)/bianchi_sim.f90 lib
$(BIANCHIBIN)/bianchi_sim:   $(BIANCHIINC)/bianchi_sim.o
	$(FC) -o $(BIANCHIBIN)/bianchi_sim $(BIANCHIINC)/bianchi_sim.o \
	$(LDFLAGS) $(PPFLAGS)



