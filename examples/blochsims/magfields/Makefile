

SHELL = /bin/sh

CC = blcc -nfmpi
CXX = blcc -nfmpi


match_SOURCES = fieldcalc.cc 
field_OBJECTS = fieldcalc.o 


INCLUDES = -I./ 

CXXCOMPILE = $(CXX) $(INCLUDES) 

SOURCES = $(rot_SOURCES)
OBJECTS = $(rot_OBJECTS)

.cc.o:
	$(CXXCOMPILE)  -c $<


all: field



field: fieldcalc.o
	$(CXXCOMPILE) $< -o fieldBF
	
	
.PHONY: 

.NOEXPORT:
