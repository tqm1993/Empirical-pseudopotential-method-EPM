
SRCDIR=trunk
MODDIR=trunk/module
MPIDIR=/usr/local/lib

FC=mpif90
FFLAGS=-ffree-form -cpp -O3  -fPIC -fno-second-underscore -m64 -Wl,--no-as-needed 

vpath %.f90 $(SRCDIR)

%.o:%.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -I$(MPIDIR)  -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lgomp -lpthread -lm -ldl -c -o $(MODDIR)/$@ $<	

.PHONY:	all	clean

all: phys_cte.o	spline.o	setup_epm.o	abs_routine.o	fit_routine.o

clean:
	@echo "Cleaning object and module files in "$(MODDIR)"."; \
	rm -f *.o $(MODDIR)/*.o	\
	rm -f *.o $(MODDIR)/*.mod
	






