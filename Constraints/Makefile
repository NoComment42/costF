SRC=class_constraints.F03
OBJ=class_constraints.o


%.o : %.f90
	$(FORT) $(FLAGS) -c $< -o $@

%.o : %.F03
	$(FORT) $(FLAGS) -c $< -o $@
	

ifeq (0,${MAKELEVEL}) 

FORT=gfortran
FLAGS= -fdefault-real-8 -O3 -fbacktrace -fcheck=all -fimplicit-none -Wall
%FLAGS= -fdefault-real-8 

SRC+= test.f90
OBJ+= test.o

all: $(OBJ)
	$(FORT) $(FLAGS) $(OBJ)

else

all: $(OBJ)

endif

clean:
	rm -f *.o *.mod

scrub: clean
	rm -f *.dat
	[ -e a.out ] || rm a.out
