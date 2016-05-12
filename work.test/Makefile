include ../make.inc

CC    = $(CPPCOMPILER) $(CPPCOMPILERFLAGS)

PROG    = testUtil testODErkck testAMod testGrad testOptim testSW testFindAtom testDispAtom testRmAtoms testRandRmAtoms testGPHPretty testGPHsVI testGPHWave testStrain2D testGPHnLa testGPHonHole testLJ testClockPot testPreClockMC testClockMC testScattHeC testMC testMD testMDRadHe testMDExx

##############################

all:    $(PROG)

testUtil: testUtil.o
	$(CC) testUtil.o -lm -o $@

testODErkck: testODErkck.o
	$(CC) testODErkck.o -lm -o $@

testAMod: testAMod.o dir_RNG dir_AMod
	$(CC) testAMod.o ../RNG/*.o ../AMod/*.o -lm -o $@

testGrad: testGrad.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII dir_donlp2
	$(CC) testGrad.o ../RNG/*.o ../AMod/*.o ../AMod/*.ko ../AMod/*.don ../REBO/*.o ../LCBOPII/*.o ../donlp2/*.o -lm -o $@

testOptim: testOptim.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII dir_donlp2
	$(CC) testOptim.o ../RNG/*.o ../AMod/*.o ../AMod/*.ko ../AMod/*.don ../REBO/*.o ../LCBOPII/*.o ../donlp2/*.o -lm -o $@

testSW: testSW.o dir_AMod dir_SW
	$(CC) testSW.o ../AMod/*.o ../SW/*.o -lm -o $@

testFindAtom: testFindAtom.o dir_AMod
	$(CC) testFindAtom.o ../AMod/*.o -lm -o $@

testDispAtom: testDispAtom.o dir_AMod
	$(CC) testDispAtom.o ../AMod/*.o -lm -o $@

testRmAtoms: testRmAtoms.o dir_AMod
	$(CC) testRmAtoms.o ../AMod/*.o -lm -o $@

testRandRmAtoms: testRandRmAtoms.o dir_AMod dir_RNG
	$(CC) testRandRmAtoms.o ../AMod/*.o ../RNG/*.o -lm -o $@

testGPHPretty: testGPHPretty.o dir_AMod
	$(CC) testGPHPretty.o ../AMod/*.o -lm -o $@

testGPHsVI: testGPHsVI.o dir_RNG dir_AMod dir_REBO dir_LCBOPII dir_donlp2
	$(CC) testGPHsVI.o ../RNG/*.o ../AMod/*.o ../AMod/*.ko ../AMod/*.don ../REBO/*.o ../LCBOPII/*.o ../donlp2/*.o -lm -o $@

testGPHWave: testGPHWave.o dir_AMod
	$(CC) testGPHWave.o ../AMod/*.o -lm -o $@

testStrain2D: testStrain2D.o dir_AMod
	$(CC) testStrain2D.o ../AMod/*.o -lm -o $@

testGPHnLa: testGPHnLa.o dir_RNG dir_AMod
	$(CC) testGPHnLa.o ../RNG/*.o ../AMod/*.o -lm -o $@

testGPHonHole: testGPHonHole.o dir_RNG dir_AMod
	$(CC) testGPHonHole.o ../RNG/*.o ../AMod/*.o -lm -o $@

testLJ: testLJ.o
	$(CC) testLJ.o -lm -o $@

testClockPot: testClockPot.o dir_AMod
	$(CC) testClockPot.o ../AMod/*.o ../AMod/*.ko ../REBO/*.o ../LCBOPII/*.o -lm -o $@

testPreClockMC: testPreClockMC.o dir_RNG
	$(CC) testPreClockMC.o ../RNG/*.o -lm -o $@

testClockMC: testClockMC.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) testClockMC.o ../RNG/*.o ../AMod/*.o ../AMod/*.ko ../MC/*.o ../REBO/*.o ../LCBOPII/*.o -lm -o $@

testScattHeC: testScattHeC.o dir_AMod dir_MC dir_LCBOPII
	$(CC) testScattHeC.o ../AMod/*.o ../MC/MC_LCBOPII.o ../LCBOPII/*.o -lm -o $@

testMC: testMC.o dir_RNG dir_AMod dir_MC dir_REBO dir_LCBOPII
	$(CC) testMC.o ../RNG/*.o ../AMod/*.o ../MC/MC_REBO.o ../REBO/*.o ../MC/MC_LCBOPII.o ../LCBOPII/*.o -lm -o $@

testMD: testMD.o  dir_RNG dir_AMod dir_REBO dir_LCBOPII dir_MD
	$(CC) testMD.o ../RNG/*.o ../AMod/*.o ../AMod/*.ko ../REBO/*.o ../LCBOPII/*.o ../MD/*.o -lm -o $@

testMDRadHe: testMDRadHe.o dir_RNG dir_AMod dir_REBO dir_LCBOPII dir_MD
	$(CC) testMDRadHe.o ../RNG/*.o ../AMod/*.o ../AMod/*.ko ../REBO/*.o ../LCBOPII/*.o ../MD/*.o -lm -o $@

testMDExx: testMDExx.o dir_RNG dir_AMod dir_REBO dir_LCBOPII dir_MD
	$(CC) testMDExx.o ../RNG/*.o ../AMod/*.o ../AMod/*.ko ../REBO/*.o ../LCBOPII/*.o ../MD/*.o -lm -o $@

cleanobj:
	$(REMOVE) *.o

clean:
	$(REMOVE) *.o $(PROG)

##############################

dir_RNG:
	( cd ../RNG && $(MAKE) )

dir_AMod:
	( cd ../AMod && $(MAKE) )

dir_SW:
	( cd ../SW && $(MAKE) )

dir_MC:
	( cd ../MC && $(MAKE) )

dir_MD:
	( cd ../MD && $(MAKE) )

dir_REBO:
	( cd ../REBO && $(MAKE) )

dir_LCBOPII:
	( cd ../LCBOPII && $(MAKE) )

dir_donlp2:
	( cd ../donlp2 && $(MAKE) )

%.o: %.cpp
	$(CC) -c -o $@ $<

