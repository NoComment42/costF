EXEC=cost
SRC=bisection.F03 cost.F03 class_constraints.F03
OBJ=bisection.o cost.o class_constraints.o
LIBDIR=FunctionClass
CONDIR=Constraints

export FORT=gfortran
export FLAGS= -cpp -fdefault-real-8 -O3 -fbacktrace -Wall -fcheck=all \
       -ffree-line-length-0
HOME= $(shell pwd)
PREPROC= -D'EXPDIR="$(HOME)/EXP"' -D'SIMDIR="."' -D'OUTDIR="$(HOME)/PRINTFILES"'
PREPROC += -D'DEBUG=1'

%.o %.mod : %.F03
	$(FORT) $(FLAGS) $(PREPROC) -c $< -o $@

	
$(EXEC): link $(OBJ)
	$(FORT) $(FLAGS) $(PREPROC) $(OBJ) -L$(LIBDIR) -lpw -o $(EXEC) 

first:
	cd $(LIBDIR) && $(MAKE)
	cd $(CONDIR) && $(MAKE)

link: first
	@[ -L class_pwfunc.mod ] \
	|| ln -s $(LIBDIR)/class_pwfunc.mod class_pwfunc.mod
	@[ -L class_pwlinear.mod ] || \
	ln -s $(LIBDIR)/class_pwlinear.mod class_pwlinear.mod
	@[ -L class_constraints.mod ] || \
	ln -s $(CONDIR)/class_constraints.mod class_constraints.mod
	@[ -L class_constraints.o ] || \
	ln -s $(CONDIR)/class_constraints.o class_constraints.o

clean:
	touch foo.o
	touch foo.mod
	rm *.o *.mod
	cd $(LIBDIR) && $(MAKE) clean
	cd $(CONDIR) && $(MAKE) clean
	
scrub:
	rm -f *.o *.mod
	rm -f $(EXEC)
	cd $(LIBDIR) && $(MAKE) scrub
	cd $(CONDIR) && $(MAKE) scrub
	
