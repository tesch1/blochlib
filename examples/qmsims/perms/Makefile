

SHELL = /bin/sh

#CC = /usr/people/magneto/packages/nompi/bin/blcc 
CXX = blcc 
#CXX =  /usr/people/magneto/packages/withmpi/bin/blcc 



permute_SOURCES = src/permute.cc src/spindex.cc
permute_OBJECTS = src/permute.o src/spindex.o

pC7_OBJECTS=src/sequenceparse.o src/pulsedata.o src/RecoupleTraingen.o src/spindex.o src/RecoupleContent.o \
	    src/RecoupleSubUnits.o
test_OBJECTS= src/testSubs.o  src/sequenceparse.o src/pulsedata.o src/RecoupleSubUnits.o
pC7_SOURCES = src/RecoupleTraingen.cc src/spindex.cc


testP_OBJECTS=src/sequenceparse.o src/pulsedata.o src/testP.o

INCLUDES = -I./ 

CXXCOMPILE = $(CXX) $(INCLUDES) 

SOURCES = $(permute_SOURCES)
OBJECTS = $(permute_OBJECTS)

.cc.o:
	$(CXXCOMPILE)  -c  $<
	mv *.o src/


all: tensorgen Recouple


tensorgen: $(permute_OBJECTS)
	$(CXXCOMPILE)  $(permute_OBJECTS) -o tensorgen


Recouple: $(pC7_OBJECTS)
	$(CXXCOMPILE)  $(pC7_OBJECTS) -o rectrain
	
test: $(test_OBJECTS)
	$(CXXCOMPILE)  $(test_OBJECTS) -o test

testP: $(testP_OBJECTS)
	$(CXXCOMPILE)  $(testP_OBJECTS) -o ptest
	

.PHONY: 

.NOEXPORT:
