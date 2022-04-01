SRCDIR   := src
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LIBS     := -lm
RM       := rm -f

EXE      := bin/mtau

LIBSRC   := $(wildcard $(SRCDIR)/*.cc)
LIBOBJ   := $(LIBSRC:.cc=.o)

# YAM2 (https://github.com/cbpark/YAM2)
YAM2     ?= /usr/local
CXXFLAGS += -I$(YAM2)/include
LIBS     += -L$(YAM2)/lib -lYAM2 -Wl,-rpath $(YAM2)/lib

# NLopt (https://nlopt.readthedocs.io/)
NLOPT    ?= /usr
LIBS     += -L$(NLOPT)/lib -lnlopt -Wl,-rpath $(NLOPT)/lib

# ROOT (https://root.cern.ch)
CXXFLAGS += -I$(shell root-config --incdir)
LIBS     += $(shell root-config --libs)
LIBS     += -lGenVector

.PHONY: all clean
.PRECIOUS: $(LIBOBJ)

all: $(EXE)

bin/%: bin/%.o $(LIBOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean::
	$(RM) $(EXE)
	$(RM) $(LIBOBJ)
