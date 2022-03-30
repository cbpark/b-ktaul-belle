SRCDIR   := src
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LIBS     := -lm
RM       := rm -f

EXE := bin/mtau

# NLopt (https://nlopt.readthedocs.io/)
NLOPT ?= /usr
LIBS  += -L$(NLOPT)/lib -lnlopt -Wl,-rpath $(NLOPT)/lib

# ROOT (https://root.cern.ch)
CXXFLAGS += -I$(shell root-config --incdir)
LIBS     += $(shell root-config --libs)

.PHONY: all clean

all: $(EXE)

bin/%: bin/%.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean::
	$(RM) $(EXE)
