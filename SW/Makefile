include ../make.inc

CC      = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

OBJS    = SWStimulus.o SWCreator.o SWRelaxor.o

##############################

all:    SW_OBJS

SW_OBJS:  $(OBJS)

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o

##############################

%.o: %.cpp
	$(CC) -c -o $@ $<

