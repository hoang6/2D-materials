include make.inc

##############################

all: dir_work.test dir_work.util dir_work.MC dir_work.phonon

dir_work.test:
	(cd work.test && $(MAKE))

dir_work.util:
	(cd work.util && $(MAKE))

dir_work.MC:
	(cd work.MC && $(MAKE))

dir_work.phonon:
	(cd work.phonon && $(MAKE))

cleanobj: 
	(cd work.test && $(MAKE) -i cleanobj)
	(cd work.util && $(MAKE) -i cleanobj)
	(cd work.MC && $(MAKE) -i cleanobj)
	(cd work.phonon && $(MAKE) -i cleanobj)
	(cd AnyOption && $(MAKE) -i cleanobj)
	(cd RNG && $(MAKE) -i cleanobj)
	(cd donlp2 && $(MAKE) -i cleanobj)
	(cd AMod && $(MAKE) -i cleanobj)
	(cd SW && $(MAKE) -i cleanobj)
	(cd MC && $(MAKE) -i cleanobj)
	(cd MD && $(MAKE) -i cleanobj)
	(cd REBO && $(MAKE) -i cleanobj)
	(cd LCBOPII && $(MAKE) -i cleanobj)

clean: 
	(cd work.test && $(MAKE) -i clean)
	(cd work.util && $(MAKE) -i clean)
	(cd work.MC && $(MAKE) -i clean)
	(cd work.phonon && $(MAKE) -i clean)
	(cd AnyOption && $(MAKE) -i clean)
	(cd RNG && $(MAKE) -i clean)
	(cd donlp2 && $(MAKE) -i clean)
	(cd AMod && $(MAKE) -i clean)
	(cd SW && $(MAKE) -i clean)
	(cd MC && $(MAKE) -i clean)
	(cd MD && $(MAKE) -i clean)
	(cd REBO && $(MAKE) -i clean)
	(cd LCBOPII && $(MAKE) -i clean)

##############################
