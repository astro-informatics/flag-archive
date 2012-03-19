# ======================================== #

# Directory for SSHT
SSHTDIR	= ${SSHT}
# Directory for FFTW
FFTWDIR	= ${FFTW}
# Directory for FFTW
GSLDIR	= ${GSL}
# Directory for MATLAB
MLAB	= /usr/local/MATLAB/R2010aSV
# Directory for DOXYGEN
DOXYGEN_PATH=/Applications/Doxygen.app/Contents/Resources/doxygen

# Compiler and options
CC	= gcc
OPT	= -Wall -O3
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

# === FLAG ===
FLAGDIR = .
FLAGLIB = $(FLAGDIR)/lib/c
FLAGINC = $(FLAGDIR)/include/c
FLAGBIN = $(FLAGDIR)/bin/c
FLAGLIBN= flag
FLAGSRC	= $(FLAGDIR)/src/c
FLAGOBJ = $(FLAGSRC)

# === SSHT ===
SSHTLIB	= $(SSHTDIR)/lib/c
SSHTINC	= $(SSHTDIR)/include/c
SSHTLIBN= ssht

# === FFTW ===
FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3

# ======================================== #

FLAGSRCMAT	= $(FLAGDIR)/src/matlab
FLAGOBJMAT  = $(FLAGSRCMAT)
FLAGOBJMEX  = $(FLAGSRCMAT)

vpath %.c $(FLAGSRC)
vpath %.h $(FLAGSRC)
vpath %_mex.c $(FLAGSRCMAT)

LDFLAGS = -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -L$(FLAGLIB) -l$(FLAGLIBN) -lm -fopenmp

LDFLAGSMEX = -I/usr/local/include -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -L$(FLAGLIB) -l$(FLAGLIBN)

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC) -I$(FLAGINC)  -fopenmp

FLAGOBJS= $(FLAGOBJ)/flag_core.o	\
	  $(FLAGOBJ)/flag_sampling.o	\
	  $(FLAGOBJ)/flag_io.o			\
	  $(FLAGOBJ)/flag_spherlaguerre.o

ifneq (,$(wildcard $(GSLDIR)/gsl_math.h))
	FLAGOBJS+= $(FLAGOBJ)/flag_spherbessel.o
	FFLAGS+= -I$(GSLDIR)
	LDFLAGS+= -L$(GSLDIR)
	LDFLAGSMEX+= -L$(GSLDIR)
endif

FLAGOBJSMAT = $(FLAGOBJMAT)/flag_analysis_mex.o	\
	$(FLAGOBJMAT)/flag_synthesis_mex.o	\
	$(FLAGOBJMAT)/flag_sampling_mex.o	\
	$(FLAGOBJMAT)/slag_synthesis_mex.o	\
	$(FLAGOBJMAT)/slag_analysis_mex.o	\
	$(FLAGOBJMAT)/slag_sampling_mex.o

FLAGOBJSMEX = $(FLAGOBJMEX)/flag_analysis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/flag_synthesis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/flag_sampling_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/slag_synthesis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/slag_analysis_mex.$(MEXEXT)	\
	$(FLAGOBJMEX)/slag_sampling_mex.$(MEXEXT)

$(FLAGOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(FLAGOBJMAT)/%_mex.o: %_mex.c $(FLAGLIB)/lib$(FLAGLIBN).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC} 

$(FLAGOBJMEX)/%_mex.$(MEXEXT): $(FLAGOBJMAT)/%_mex.o $(FLAGLIB)/lib$(FLAGLIBN).a
	$(MEX) $< -o $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

# ======================================== #

.PHONY: default
default: lib test tidy

.PHONY: matlab
matlab: $(FLAGOBJSMEX)

.PHONY: all
all: lib test matlab tidy

.PHONY: lib
lib: $(FLAGLIB)/lib$(FLAGLIBN).a
$(FLAGLIB)/lib$(FLAGLIBN).a: $(FLAGOBJS)
	ar -r $(FLAGLIB)/lib$(FLAGLIBN).a $(FLAGOBJS)

.PHONY: test
test: $(FLAGBIN)/flag_test
$(FLAGBIN)/flag_test: $(FLAGOBJ)/flag_test.o $(FLAGLIB)/lib$(FLAGLIBN).a
	$(CC) $(OPT) $< -o $(FLAGBIN)/flag_test $(LDFLAGS)

.PHONY: doc
doc:
	$(DOXYGEN_PATH) $(FLAGDIR)/src/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(FLAGDIR)/doc/html/*

.PHONY: clean
clean:	tidy
	rm -f $(FLAGLIB)/lib$(FLAGLIBN).a
	rm -f $(FLAGOBJMEX)/*_mex.$(MEXEXT)
	rm -f $(FLAGBIN)/flag_test

.PHONY: tidy
tidy:
	rm -f $(FLAGOBJ)/*.o
	rm -f *~ 

