CLEANFILES = WJetsGenerator.log WJetsGenerator.out \
             WJetsGenerator.run WJetsGenerator.tex EETest.log EETest.out \
             EETest.run EETest.tex done-all-links \
             W0cache.dat.gz W1cache.dat.gz W2cache.dat.gz W3cache.dat.gz \
             .runThePEG.timer.WJetsGenerator.run LHETest.log LHETest.out \
             LHETest.run LHETest.tex LHETest1.log LHETest1.out \
             LHETest1.run LHETest1.tex LHETest2.log LHETest2.out \
             LHETest2.run LHETest2.tex LHETest3.log LHETest3.out \
             LHETest3.run LHETest3.tex LHETest4.log LHETest4.out \
             LHETest4.run LHETest4.tex LHETest5.log LHETest5.out \
             LHETest5.run LHETest5.tex .runThePEG.timer.LHETest.run \
             W1jetsCache.gz W1jetsuCache.gz W1jetsMaxCache.gz \
             W1jetsuMaxCache.gz W1jetsMinCache.gz W1jetsuMinCache.gz \
             W0jetsCache.gz W0jetsuCache.gz W0jetsMaxCache.gz \
             W0jetsuMaxCache.gz W0jetsMinCache.gz W0jetsuMinCache.gz

EXTRA_DIST = WJetsGenerator.in W0jets.dat.gz W1jets.dat.gz W2jets.dat.gz \
             W3jets.dat.gz LHETest.in LHCW0J.dat.gz LHCW0Ju.dat.gz \
             LHCW1J.dat.gz LHCW1Ju.dat.gz EETest.in ee4jet.lhe.gz

INPUTFILES = WJetsGenerator.in W0jets.dat.gz W1jets.dat.gz W2jets.dat.gz \
             W3jets.dat.gz LHETest.in LHCW0J.dat.gz LHCW0Ju.dat.gz \
             LHCW1J.dat.gz LHCW1Ju.dat.gz EETest.in ee4jet.lhe.gz

all-local: done-all-links

done-all-links:
@EMPTY@ifdef SHOWCOMMAND
	for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	echo "stamp" > done-all-links
@EMPTY@else
	@echo "sym-linking input files files..."
	@for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	@echo "stamp" > done-all-links
@EMPTY@endif

WJetsGenerator.run: done-all-links WJetsGenerator.in W0jets.dat.gz W1jets.dat.gz W2jets.dat.gz W3jets.dat.gz
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo WJetsGenerator.in 

WJetsGenerator.out: WJetsGenerator.run
	time $(RUNTHEPEG) -d 1 WJetsGenerator.run

EETest.run: done-all-links EETest.in ee4jet.lhe.gz
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo EETest.in 

EETest.out: EETest.run
	time $(RUNTHEPEG) -d 1 EETest.run

if USELHAPDF
LHETest.run: done-all-links LHETest.in LHCW0J.dat.gz LHCW0Ju.dat.gz \
             LHCW1J.dat.gz LHCW1Ju.dat.gz
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo LHETest.in 

LHETest.out: LHETest.run
	time $(RUNTHEPEG) -d 1 LHETest.run
endif

check-local:
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo WJetsGenerator.in 
	time $(RUNTHEPEG) -d 1 WJetsGenerator.run
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo EETest.in 
	time $(RUNTHEPEG) -d 1 EETest.run
if USELHAPDF
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo LHETest.in 
	time $(RUNTHEPEG) -d 1 LHETest.run
endif

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(LN_S) ../ee4jet.lhe.gz ../W0jets.dat.gz ../W1jets.dat.gz ../W2jets.dat.gz ../W3jets.dat.gz .; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../WJetsGenerator.in; $(DESTDIR)$(bindir)/runThePEG WJetsGenerator.run; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../EETest.in; $(DESTDIR)$(bindir)/runThePEG EETest.run
if USELHAPDF
	cd _test; $(LN_S)  ../LHCW0J.dat.gz ../LHCW0Ju.dat.gz ../LHCW1J.dat.gz ../LHCW1Ju.dat.gz .; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../LHETest.in; $(DESTDIR)$(bindir)/runThePEG LHETest.run
endif
	rm -rf _test

save:
	mkdir -p save
	cp *.log *.out *.run *.tex save


.PHONY: save
