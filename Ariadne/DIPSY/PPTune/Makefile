%: 
	./expand.pl -p $@ PPTune.inc > PPTune.in
	../../../ThePEG/src/setupThePEG -L ../../../ThePEG/lib -L../../../Pythia7/lib -L../../../TheP8I/lib -L../../lib --exitonerror -L ../.libs -r ../DIPSYDefaults.rpo PPTune.in

