# ======================================== #

# Directory for SSHT
SSHTDIR	= ${SSHT}
# Directory for FFTW
FFTWDIR	= ${FFTW}
# Directory for GSL
GSLDIR	= ${GSL}
# Directory for MATLAB
MLAB	=  /Applications/MATLAB_R2011b.app
# Directory for DOXYGEN
DOXYGEN_PATH=doxygen

# Compiler and options
CC	= gcc
OPT	= -Wall -O3 -g -DFLAG_VERSION=\"1.0\" -DFLAG_BUILD=\"`svnversion -n .`\"
UNAME := $(shell uname)

# ======================================== #

# === MATLAB ===
ifeq ($(UNAME), Linux)
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib
  MEXEXT	= mexa64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif
ifeq ($(UNAME), Darwin)
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib
  MEXEXT	= mexmaci64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif

# === SSHT ===
GSLLIB	= $(GSLDIR)/lib
GSLINC	= $(GSLDIR)/include/gsl
GSLLIBN= gsl

# === FLAG ===
FLAGDIR = .
FLAGLIB = $(FLAGDIR)/lib
FLAGINC = $(FLAGDIR)/include
FLAGBIN = $(FLAGDIR)/bin
FLAGLIBN= flag
FLAGSRC	= $(FLAGDIR)/src/main/c
FLAGTESTSRC	= $(FLAGDIR)/src/test/c
FLAGOBJ = $(FLAGSRC)
FLAGTESTOBJ = $(FLAGTESTSRC)

# === SSHT ===
SSHTLIB	= $(SSHTDIR)/lib/c
SSHTINC	= $(SSHTDIR)/include/c
SSHTLIBN= ssht

# === FFTW ===
FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3

# ======================================== #

FLAGSRCMAT	= $(FLAGDIR)/src/main/matlab
FLAGOBJMAT  = $(FLAGSRCMAT)
FLAGOBJMEX  = $(FLAGSRCMAT)

vpath %.c $(FLAGSRC)
vpath %.c $(FLAGTESTSRC)
vpath %.h $(FLAGINC)
vpath %_mex.c $(FLAGSRCMAT)

LDFLAGS = -L$(FLAGLIB) -l$(FLAGLIBN) -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -lm -lc -fopenmp

LDFLAGSMEX = -L$(FLAGLIB) -l$(FLAGLIBN) -I/usr/local/include -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN)

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC) -I$(FLAGINC)  -fopenmp

FLAGOBJS= $(FLAGOBJ)/flag_core.o	\
	  $(FLAGOBJ)/flag_sampling.o	\
	  $(FLAGOBJ)/flag_io.o			\
	  $(FLAGOBJ)/flag_spherlaguerre.o			\
	  $(FLAGOBJ)/flag_spherbessel.o

ifneq (,$(wildcard $(GSLINC)/gsl_sf.h))
	FLAGOBJS+= $(FLAGOBJ)/flag_spherbessel.o
	FFLAGS+= -I$(GSLINC)
	LDFLAGS+= -L$(GSLLIB)
	LDFLAGS+= -l$(GSLLIBN)
	LDFLAGSMEX+= -L$(GSLLIB)
	LDFLAGSMEX+= -l$(GSLLIBN)
endif

FLAGOBJSMAT = $(FLAGOBJMAT)/flag_analysis_mex.o	\
	$(FLAGOBJMAT)/flag_synthesis_mex.o	\
	$(FLAGOBJMAT)/flag_sampling_mex.o	\
	$(FLAGOBJMAT)/slag_synthesis_mex.o	\
	$(FLAGOBJMAT)/slag_analysis_mex.o	\
	$(FLAGOBJMAT)/slag_sampling_mex.o	\
	$(FLAGOBJMAT)/flag_get_tau_mex.o	\
	$(FLAGOBJMAT)/slag_gausslaguerre_quadrature_mex.o

FLAGOBJSMEX = $(FLAGOBJMEX)/flag_analysis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/flag_synthesis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/flag_sampling_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/slag_synthesis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/slag_analysis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/slag_sampling_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/flag_get_tau_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/slag_gausslaguerre_quadrature_mex.$(MEXEXT)

$(FLAGOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(FLAGTESTOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(FLAGOBJMAT)/%_mex.o: %_mex.c $(FLAGLIB)/lib$(FLAGLIBN).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC} 

$(FLAGOBJMEX)/%_mex.$(MEXEXT): $(FLAGOBJMAT)/%_mex.o $(FLAGLIB)/lib$(FLAGLIBN).a
	$(MEX) $< -o $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

# ======================================== #

.PHONY: default
default: lib test about tidy

.PHONY: matlab
matlab: lib $(FLAGOBJSMEX) about

.PHONY: all
all: lib matlab doc test about tidy

.PHONY: lib
lib: $(FLAGLIB)/lib$(FLAGLIBN).a
$(FLAGLIB)/lib$(FLAGLIBN).a: $(FLAGOBJS)
	ar -r $(FLAGLIB)/lib$(FLAGLIBN).a $(FLAGOBJS)

.PHONY: test
test: lib $(FLAGBIN)/flag_test
$(FLAGBIN)/flag_test: $(FLAGTESTOBJ)/flag_test.o $(FLAGLIB)/lib$(FLAGLIBN).a
	$(CC) $(OPT) $< -o $(FLAGBIN)/flag_test $(LDFLAGS)

.PHONY: fbtest
fbtest: lib $(FLAGBIN)/flag_fbtest
$(FLAGBIN)/flag_fbtest: $(FLAGTESTOBJ)/flag_fbtest.o $(FLAGLIB)/lib$(FLAGLIBN).a
	$(CC) $(OPT) $< -o $(FLAGBIN)/flag_fbtest $(LDFLAGS)

.PHONY: about
about: $(FLAGBIN)/flag_about
$(FLAGBIN)/flag_about: $(FLAGOBJ)/flag_about.o 
	$(CC) $(OPT) $< -o $(FLAGBIN)/flag_about
	$(FLAGBIN)/flag_about

.PHONY: doc
doc:
	$(DOXYGEN_PATH) $(FLAGDIR)/src/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(FLAGDIR)/doc/c/*

.PHONY: clean
clean:	tidy cleandoc
	rm -f $(FLAGLIB)/lib$(FLAGLIBN).a
	rm -f $(FLAGOBJMEX)/*_mex.$(MEXEXT)
	rm -f $(FLAGBIN)/flag_fbtest
	rm -f $(FLAGBIN)/flag_test
	rm -f $(FLAGBIN)/flag_about

.PHONY: tidy
tidy:
	rm -f $(FLAGOBJ)/*.o
	rm -f $(FLAGOBJMEX)/*.o
	rm -f *~ 

# ======================================== #