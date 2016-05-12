include ../make.inc

CC      = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

PROG    = runMCNVT runMCNPT runMCNPH runMCNtH analMCNVT analMCNPH analMCNtH analGPH analCNT

##############################

all:    $(PROG)

runMCNVT: runMCNVT.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) runMCNVT.o ../RNG/*.o ../AMod/*.o ../MC/MC*.o ../MC/Stream.o ../MC/alloc.o ../MC/useNVT.o ../REBO/*.o ../LCBOPII/*.o $(LIB) -lm -o $@

runMCNPT: runMCNPT.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) runMCNPT.o ../RNG/*.o ../AMod/*.o ../MC/MC*.o ../MC/Stream.o ../MC/alloc.o ../MC/useNPT.o ../REBO/*.o ../LCBOPII/*.o $(LIB) -lm -o $@

runMCNPH: runMCNPH.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) runMCNPH.o ../RNG/*.o ../AMod/*.o ../MC/MC*.o ../MC/Stream.o ../MC/alloc.o ../MC/useNPH.o ../REBO/*.o ../LCBOPII/*.o $(LIB) -lm -o $@

runMCNtH: runMCNtH.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) runMCNtH.o ../RNG/*.o ../AMod/*.o ../MC/MC*.o ../MC/Stream.o ../MC/alloc.o ../MC/useNtH.o ../REBO/*.o ../LCBOPII/*.o $(LIB) -lm -o $@

analMCNVT: analMCNVT.o
	$(CC) analMCNVT.o -lm -o $@

analMCNPH: analMCNPH.o
	$(CC) analMCNPH.o -lm -o $@

analMCNtH: analMCNtH.o
	$(CC) analMCNtH.o -lm -o $@

analGPH: analGPH.o dir_AMod
	$(CC) analGPH.o ../AMod/*.o -lm -o $@

analCNT: analCNT.o dir_AMod
	$(CC) analCNT.o ../AMod/*.o -lm -o $@

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o $(PROG)

##############################

dir_RNG:
	( cd ../RNG && $(MAKE) )

dir_AMod:
	( cd ../AMod && $(MAKE) )

dir_MC:
	( cd ../MC && $(MAKE) )

dir_REBO:
	( cd ../REBO && $(MAKE) )

dir_LCBOPII:
	( cd ../LCBOPII && $(MAKE) )

%.o: %.cpp
	$(CC) -c -o $@ $<
