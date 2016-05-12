include ../make.inc

CC      = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

MPICC   = $(MPICOMPILER) $(MPICOMPILERFLAGS) 

PROG = nRecords sampleRecords manufacture 
ifneq ($(MPICOMPILER),)
	MPIPROG = simpleMPIDriver
endif

##############################

all:    $(PROG) $(MPIPROG)

nRecords: nRecords.o
	$(CC) nRecords.o -lm -o $@

sampleRecords: sampleRecords.o
	$(CC) sampleRecords.o -lm -o $@

manufacture: manufacture.o dir_AnyOption dir_AMod
	$(CC) manufacture.o ../AnyOption/*.o ../AMod/*.o -lm -o $@

simpleMPIDriver: simpleMPIDriver.mpio
	$(MPICC) simpleMPIDriver.mpio -lm -o $@

cleanobj:
	$(REMOVE) *.o *.mpio

clean:
	$(REMOVE) *.o *.mpio $(PROG) $(MPIPROG)

##############################

dir_AnyOption:
	(cd ../AnyOption && $(MAKE) )

dir_RNG:
	( cd ../RNG && $(MAKE) )

dir_AMod:
	( cd ../AMod && $(MAKE) )

%.o: %.cpp
	$(CC) -c -o $@ $<

%.mpio: %.cpp
	$(MPICC) -c -o $@ $<	
