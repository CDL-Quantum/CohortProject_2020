.SUFFIXES: .o .c .f .cxx
# DEBUG option remove for efficiency
CC=g++
cc=gcc
FC=gfortran
options= -DDUSHQUICK -DBLAS -DLINUX -O2 -Wno-deprecated
objects=FC.o matvec.o prod.o
dependencies=matvec.h 
FC_fast : $(objects)
	$(CC) -o FCF_calc $(objects)  \
-lblas -llapack -lgfortran
$(objects): $(dependencies)
.cxx.o:
	$(CC) -c $(options) $<
.c.o:
	$(cc) -c -O2 $<
.f.o:
	$(FC) -c -O2 $<
