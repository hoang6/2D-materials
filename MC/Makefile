include ../make.inc

CC    = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

OBJS    = MC_REBO.o MC_LCBOPII.o MCAgent.o MCAgent_CtrlZ.o MCAgent_STube.o Stream.o useNVT.o useNPT.o useNPH.o useNtH.o alloc.o

##############################

all: MC_OBJS

MC_OBJS: $(OBJS)

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o

##############################

%.o: %.cpp
	$(CC) -c -o $@ $<

