

SHELL = /bin/sh

CC = blcc -nmpi
CXX = blcc -nmpi


best_SOURCES = best_bfield.cc Dcoil.cc magfitter.cc exsituSim.cc
best_OBJECTS = best_bfield.o Dcoil.o magfitter.o exsituSim.o
rot_SOURCES =Dcoilcalc.cc Dcoil.cc
rot_OBJECTS =Dcoilcalc.o  Dcoil.o

$LD_LIBRARY_PATH=/usr/people/magneto/lib
INCLUDES = -I./ 

CXXCOMPILE = $(CXX) $(INCLUDES) 

SOURCES = $(rot_SOURCES)
OBJECTS = $(rot_OBJECTS)

.cc.o:
	$(CXXCOMPILE)  -c $<


all: fitter


fitter: $(best_OBJECTS)
	$(CXXCOMPILE) $(best_OBJECTS) -Wl,--rpath -Wl,/usr/people/magneto/lib  -o matchBF
	
field: $(rot_OBJECTS)
	$(CXXCOMPILE)  $(rot_OBJECTS) -Wl,--rpath -Wl,/usr/people/magneto/lib  -o fieldBF
	
.PHONY: 

.NOEXPORT:
