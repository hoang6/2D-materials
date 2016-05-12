include ../make.inc

# use the compiler call appropriate for your platform
CC = $(CCOMPILER) $(CCOMPILERFLAGS)
# note that ANSI-C is required to compile donlp2.c you may have
# to add appropriate options in CFLAGS for your native cc or you may use
# gcc


# you may add or alter compiler options  below -O3 means full optimization
#due to stacksize limitations for large n and n+nlin+nonlin  you might be forced
# to change o8para.h

BASIS = *.h

# for using the old fashioned function form take 
# DRIVER = userfu.o user_eval.o user_init.o  with the codes from
# EXAMPLES_OLD

DRIVER = user_eval.o

OPTIMIZER = donlp2.o 

OBJS = $(DRIVER) $(OPTIMIZER)

##############################

all: DONLP2_OBJS

DONLP2_OBJS: $(OBJS)

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o

##############################

.c.o:
	$(CC) $(CFLAGS) -c $<
donlp2.o : donlp2.c $(BASIS)
user_eval.o : user_eval.c $(BASIS)

