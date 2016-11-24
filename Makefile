# Use version 4.9 of g++
CXX=g++
CC=$(CXX)

# Set this to include SeqAn libraries, either system wide
# or download into current folder and set to .
SEQAN_LIB=.

CXXFLAGS+=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1 -std=c++14
LDLIBS=-lz -lpthread
DATE=on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
CXXFLAGS+=-DDATE=\""$(DATE)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result

# DEBUG build
#CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
#CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0

all: CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
all: BAMQC

debug: CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1
debug: BAMQC

BAMQC:BAMQC.o

BAMQC.o: BAMQC.cpp BAMQC.h parse.h

clean:
	rm -f *.o BAMQC
