CC = gfortran
LDFLAGS = -L/opt/local/lib
LIBS = -lnetcdff -llapack
OMP = -fopenmp

OBJS = read_qgpv.o
OMPOBJS = sor_m.o
FILES = inv_qgpv

.SUFFIXES: .f90 .o

.f90.o:
	$(CC) $(LDFLAGS) $(LIBS) -c $<

all:	$(FILES)

inv_qgpv:	$(OBJS) $(OMPOBJS) inv_qgpv.o
	$(CC) $(LDFLAGS) $(LIBS) $(OMP) -o $@ $^

$(OBJS):	read_qgpv.f90 
$(OMPOBJS):	sor_m.f90
	$(CC) $(LDFLAGS) $(LIBS) $(OMP) -c $<
inv_qgpv.o:	inv_qgpv.f90

clean:;		rm -f *.o *.mod inv_qgpv psi.grd