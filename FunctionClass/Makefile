SRC=class_pwFunc.F03 class_pwLinear.F03
OBS=class_pwFunc.o class_pwLinear.o

%.o %.mod : %.f03
	$(FORT) $(FLAGS) -c $< -o $@

%.o %.mod : %.F03
	$(FORT) $(FLAGS) -c $< -o $@
	

ifeq (0,${MAKELEVEL})

FORT=gfortran
FLAGS= -fdefault-real-8

SRC+= test.f03
OBS+= test.o

all : $(OBS)
	ar cr libpw.a *.o
	$(FORT) $(FLAGS) test.o -L./ -lpw

else

all : $(OBS)
	ar cr libpw.a *.o

endif

clean:
	rm -f *.o *.mod
	rm -f libpw.a
	
scrub: clean
	[ -e a.out ] || rm a.out
	
