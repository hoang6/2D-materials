include ../make.inc

CC    = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

OBJS    = Data.o NVE.o NVT.o NVE1He.o IO.o Sim.o

##############################

all: MD_OBJS

MD_OBJS: $(OBJS)

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o

##############################

%.o: %.cpp
	$(CC) -c -o $@ $<

