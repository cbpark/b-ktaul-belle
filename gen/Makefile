CXXFLAGS := -g -O2 -Wall -Wextra -std=c++17 -pedantic $(CXXFLAGS)
LIBS     :=
RM       := rm -f

EXE      := gen_sig_bb

# Pythia 8
CXXFLAGS += $(shell pythia8-config --cxxflags)
LIBS     += $(shell pythia8-config --libs)

# ROOT
CXXFLAGS += -I$(shell root-config --incdir)
LIBS     += $(shell root-config --libs)

.PHONY: all clean

all: $(EXE)

gen_sig_bb: gen_sig_bb.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean::
	$(RM) gen_sig_bb gen_sig_bb.o
