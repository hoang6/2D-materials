include ../make.inc

CC    = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

OBJS = AnyOption.o ParseMainArgs.o

##############################

all:    AnyOption_OBJS

AnyOption_OBJS: $(OBJS)

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o

##############################

%.o: %.cpp
	$(CC) -c -o $@ $<
