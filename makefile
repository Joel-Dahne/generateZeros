# The location of your CXSC directory
# THIS NEEDS TO BE CHANGED TO MATCH YOUR SYSTEM
CXSCDIR=/home/urathai/cxsc

LINK_TARGETS = build/generateZeros build/integrate
OBJS = build/generateZeros.o build/integrate.o build/citaylor.o	\
build/function.o
REBUILDABLES = $(OBJS) $(LINK_TARGETS)

# Optional flags to give to the C++ compiler
CXXOPTS=-Wall -Winline

# Additional include path
CXSCINC=-I$(CXSCDIR)/include -L$(CXSCDIR)/lib

# Flags to give to the C++ compiler
CXSCFLAGS=$(CXSCINC) $(CXXOPTS)

RPATH=-Wl,-R$(CXSCDIR)/lib

clean:
	rm -f $(REBUILDABLES)

all: $(LINK_TARGETS)

build/generateZeros: build/generateZeros.o build/citaylor.o build/function.o
	g++ -o $@ $(CXSCINC) $(RPATH) $^ -lcxsc

build/integrate: build/integrate.o build/citaylor.o build/function.o
	g++ -o $@ -fopenmp $(CXSCINC) $(RPATH) $^ -lcxsc

build/%.o: src/%.cpp | build
	g++ -fopenmp -o $@ $(CXSCINC) $(RPATH) -c $< -lcxsc

build:
	mkdir build
