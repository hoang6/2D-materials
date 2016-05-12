include ../make.inc

CC      = $(CCOMPILER) $(CCOMPILERFLAGS) $(RNGFLAGS)

OBJS    = SFMT.o

##############################	

all:    RNG_OBJS

RNG_OBJS:    $(OBJS)

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o

##############################

SFMT.o: SFMT.c
	$(CC) -c -o $@ $<
