include ../make.inc

CC      = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

OBJS    = Atom.o Box.o Molecule.o MolTopo.o MolIO.o MolFactory.o MolOptim.don MolAdjuster.o MolScaler.o MolAnal.o MolChecker.o PTE.o Trans.o KernelREBO.ko KernelLCBOPIIN.ko

##############################

all:    AMod_OBJS

AMod_OBJS:  $(OBJS)

cleanobj:
	$(REMOVE) *.o *.don *.ko

clean:
	$(REMOVE) *.o *.don *.ko

##############################

%.o: %.cpp
	$(CC) -c -o $@ $<

MolOptim.don: MolOptim.cpp
	$(CC) -c -o $@ $<

KernelREBO.ko: KernelREBO.cpp
	$(CC) -c -o $@ $<

KernelLCBOPIIN.ko: KernelLCBOPIIN.cpp
	$(CC) -c -o $@ $<;
