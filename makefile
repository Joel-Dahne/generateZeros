# The location of your CXSC directory
# THIS NEEDS TO BE CHANGED TO MATCH YOUR SYSTEM
CXSCDIR=/home/urathai/cxsc

LINK_TARGETS = build/generateZeros build/integrate
OBJS = build/generateZeros.o build/integrate.o build/citaylor.o	\
build/function.o
REBUILDABLES = $(OBJS) $(LINK_TARGETS)

# Optional flags to give to the C++ compiler
CXXOPTS=-O3 -fopenmp

# Additional include path
CXSCINC=-I$(CXSCDIR)/include -L$(CXSCDIR)/lib

# Flags to give to the C++ compiler
CXSCFLAGS=$(CXXOPTS) $(CXSCINC)

RPATH=-Wl,-R$(CXSCDIR)/lib

all: $(LINK_TARGETS)

clean:
	rm -f $(REBUILDABLES)

build/generateZeros: build/generateZeros.o build/citaylor.o build/function.o
	g++ -o $@ $(CXSCFLAGS) $(RPATH) $^ -lcxsc

build/integrate: build/integrate.o build/citaylor.o build/function.o
	g++ -o $@ $(CXSCFLAGS) $(RPATH) $^ -lcxsc

build/%.o: src/%.cpp | build
	g++ -o $@ $(CXSCFLAGS) $(RPATH) -c $< -lcxsc

build:
	mkdir build
