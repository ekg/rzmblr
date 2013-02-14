CXX=g++
BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib
FASTAHACK = fastahack/Fasta.o fastahack/split.o
SMITHWATERMAN = smithwaterman/SmithWatermanGotoh.o
SMITHWATERMANOBS = smithwaterman/SmithWatermanGotoh.o \
	smithwaterman/disorder.c \
	smithwaterman/IndelAllele.o \
	smithwaterman/LeftAlign.o \
	smithwaterman/Repeats.o
CXXFLAGS=-O3 -I$(BAMTOOLS_ROOT)/include -L./ #-L$(BAMTOOLS_ROOT)/lib

all: rzmblr

# builds bamtools static lib, and copies into root
libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && make
	cp bamtools/lib/libbamtools.a ./

# statically compiles bamaddrg against bamtools static lib
rzmblr: rzmblr.cpp libbamtools.a Read.o ssw_cpp.o ssw.o $(FASTAHACK) $(SMITHWATERMAN)
	$(CXX) $(CXXFLAGS) rzmblr.cpp Read.o ssw_cpp.o ssw.o $(FASTAHACK) $(SMITHWATERMANOBS) -Ismithwaterman/ -Ibamtools/src/ -o rzmblr -lbamtools -lz

Read.o: Read.h Read.cpp
ssw.o:ssw.h
ssw_cpp.o:ssw_cpp.h

$(FASTAHACK):
	cd fastahack && $(MAKE)

$(SMITHWATERMAN):
	cd smithwaterman && $(MAKE)

clean:
	rm -f rzmblr libbamtools.a *.o
	cd fastahack && $(MAKE) clean
	cd smithwaterman && $(MAKE) clean

.PHONY: clean
