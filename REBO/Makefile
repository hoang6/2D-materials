include ../make.inc

CC      = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

OBJS    = Energy.o Force.o constants.o theta2.o g.o pij.o radic.o

##############################

all: REBO_OBJS

REBO_OBJS: $(OBJS)

cleanobj:
	$(REMOVE) *.o inHeaderMaker/in2HeaderMaker inHeaderMaker/in3HeaderMaker

clean:
	$(REMOVE) *.o inHeaderMaker/in2HeaderMaker inHeaderMaker/in3HeaderMaker

##############################

%.o: %.cpp
	$(CC) -c -o $@ $<

radic.o: in3.h

in3.h: inHeaderMaker/in3HeaderMaker
	cd inHeaderMaker && $(CURRENTDIR)in3HeaderMaker;

inHeaderMaker/in3HeaderMaker: inHeaderMaker/in3HeaderMaker.cpp
	cd inHeaderMaker && $(CC) -o in3HeaderMaker in3HeaderMaker.cpp 

pij.o: in2.h

in2.h: inHeaderMaker/in2HeaderMaker
	cd inHeaderMaker && $(CURRENTDIR)in2HeaderMaker;

inHeaderMaker/in2HeaderMaker: inHeaderMaker/in2HeaderMaker.cpp
	cd inHeaderMaker && $(CC) -o in2HeaderMaker in2HeaderMaker.cpp 

