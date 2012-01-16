
CC	= gcc
OPT	= -Wall -O3

UNAME := $(shell uname)

FLAGDIR = .
FLAGLIB = $(FLAGDIR)/lib/c
FLAGINC = $(FLAGDIR)/include/c
FLAGBIN = $(FLAGDIR)/bin/c
FLAGLIBN= flag
FLAGSRC	= $(FLAGDIR)/src/c
FLAGOBJ = $(FLAGSRC)

SSHTDIR	= ${SSHT}
SSHTLIB	= $(SSHTDIR)/lib/c
SSHTINC	= $(SSHTDIR)/include/c
SSHTLIBN= ssht

vpath %.c $(FLAGSRC)
vpath %.h $(FLAGSRC)

LDFLAGS = -L$(SSHTLIB) -l$(SSHTLIBN) -L$(FLAGLIB) -l$(FLAGLIBN) -lm

FFLAGS  = -I$(SSHTINC) -I$(FLAGINC) 

FLAGOBJS= $(FLAGOBJ)/flag_core.o	\
	  $(FLAGOBJ)/flag_sampling.o	\
	  $(FLAGOBJ)/flag_spherbessel.o	\
	  $(FLAGOBJ)/flag_spherlaguerre.o

FLAGHEADERS = flag_spherlaguerre.h

$(FLAGOBJ)/%.o: %.c $(FLAGHEADERS)
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

.PHONY: default
default: lib test tidy

.PHONY: about
about: $(FLAGBIN)/ssht_about
$(FLAGBIN)/flag_about: $(FLAGOBJ)/flag_about.o 
	$(CC) $(OPT) $< -o $(FLAGBIN)/flag_about

.PHONY: lib
lib: $(FLAGLIB)/lib$(FLAGLIBN).a
$(FLAGLIB)/lib$(FLAGLIBN).a: $(FLAGOBJS)
	ar -r $(FLAGLIB)/lib$(FLAGLIBN).a $(FLAGOBJS)

.PHONY: test
lib: $(FLAGBIN)/flag_test
$(FLAGBIN)/flag_test: $(FLAGOBJ)/flag_test.o $(FLAGLIB)/lib$(FLAGLIBN).a
	$(CC) $(OPT) $< -o $(FLAGBIN)/flag_test $(LDFLAGS)

.PHONY: clean
clean:	tidy
	rm -f $(FLAGOBJ)/*.o
	rm -f $(FLAGLIB)/lib$(FLAGLIBN).a

.PHONY: tidy
tidy:
	rm -f *~ 

