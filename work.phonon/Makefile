include ../make.inc

CC      = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

PROG    = makeGPHGroundState getGPHDynMat getGPHFreeEnergy

##############################

all:    $(PROG)

makeGPHGroundState: makeGPHGroundState.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) makeGPHGroundState.o ../RNG/*.o ../AMod/*.o ../MC/Stream.o ../MC/MC_REBO.o ../MC/MC_LCBOPII.o ../REBO/*.o ../LCBOPII/*.o $(LIB) -lm -o $@

getGPHDynMat: getGPHDynMat.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) getGPHDynMat.o ../RNG/*.o ../AMod/*.o ../MC/Stream.o ../MC/MC_REBO.o ../MC/MC_LCBOPII.o ../REBO/*.o ../LCBOPII/*.o $(LIB) -lm -o $@

getGPHFreeEnergy: getGPHFreeEnergy.o dir_AMod
	$(CC) getGPHFreeEnergy.o ../AMod/*.o $(LIB) -lm -o $@

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
