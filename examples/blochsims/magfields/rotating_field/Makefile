

SHELL = /bin/sh

CC = blcc -nfmpi
CXX = blcc -nfmpi


match_SOURCES = match_bfield.cc 
field_OBJECTS = match_bfield.o 

rot_SOURCES =blochsim_rot_bfield.cc blochsim_step_bfield.cc
rot_OBJECTS =blochsim_rot_bfield.o blochsim_step_bfield.o


INCLUDES = -I./ 

CXXCOMPILE = $(CXX) $(INCLUDES) 

SOURCES = $(rot_SOURCES)
OBJECTS = $(rot_OBJECTS)

.cc.o:
	$(CXXCOMPILE)  -c $<


all: match rot1 rot2


match: match_bfield.o
	$(CXXCOMPILE) $< -o matchBF
	
rot1: blochsim_rot_bfield.o
	$(CXXCOMPILE)  $< -o rotateBF

rot2: blochsim_step_bfield.o
	$(CXXCOMPILE)  $< -o stepBF
	
.PHONY: 

.NOEXPORT:
