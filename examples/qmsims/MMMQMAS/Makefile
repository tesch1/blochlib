

SHELL = /bin/sh

CC = blcc -n
CXX = blcc -n


reduce_SOURCES = reduce.cc propreduce.cc
reduce_OBJECTS = reduce.o propreduce.o

mmqmas_SOURCES =propreduce.cc mmmqmas.cc
mmqmas_OBJECTS =propreduce.o mmmqmas.o


INCLUDES = -I./ 

CXXCOMPILE = $(CXX) $(INCLUDES) 

SOURCES = $(mmqmas_SOURCES)
OBJECTS = $(mmqmas_OBJECTS)

.cc.o:
	$(CXXCOMPILE)  -c $<


all: reduce mmqmas


reduce: $(reduce_OBJECTS)
	$(CXXCOMPILE)  $(reduce_OBJECTS) -o reduce


mmqmas: $(mmqmas_OBJECTS)
	$(CXXCOMPILE)  $(mmqmas_OBJECTS) -o mmqmas
.PHONY: 

.NOEXPORT:
