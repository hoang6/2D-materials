include ../make.inc

CC      = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

OBJS    = Energy.o b.o FAT.o SVsr.o SVmr.o SVlr.o

##############################

all: LCBOPII_OBJS

LCBOPII_OBJS: $(OBJS)

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o

##############################

%.o: %.cpp
	$(CC) -c -o $@ $<
