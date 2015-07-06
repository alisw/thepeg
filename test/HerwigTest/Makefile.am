AUTOMAKE_OPTIONS = -Wno-portability

CLEANFILES = TestHerwig.log TestHerwig.out TestHerwig.run TestHerwig.tex \
             TestHerwig.aida .runThePEG.timer.TestHerwig.run \
             MadTest.log MadTest.out MadTest.run MadTest.tex \
             .runThePEG.timer.MadTest.run \
             done-all-links

pkglib_LTLIBRARIES = MEqq2gZ2ll.la DYMass.la

# Version info should be updated if any interface or persistent I/O
# function is changed
MEqq2gZ2ll_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
MEqq2gZ2ll_la_SOURCES = MEqq2gZ2ll.cc MEqq2gZ2ll.h MEqq2gZ2ll.icc

# Version info should be updated if any interface or persistent I/O
# function is changed
DYMass_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
DYMass_la_SOURCES = DYMass.cc DYMass.h DYMass.fh DYMass.icc

EXTRA_DIST = TestHerwig.in MadTest.in pp2mumu.dat.gz

INPUTFILES = TestHerwig.in MadTest.in pp2mumu.dat.gz

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

TestHerwig.run: done-all-links TestHerwig.in 
	$(SETUPTHEPEG) --exitonerror -L .libs -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo TestHerwig.in 


TestHerwig.out: TestHerwig.run MEqq2gZ2ll.la DYMass.la
	time $(RUNTHEPEG) -L .libs -d 1 TestHerwig.run

MadTest.run: done-all-links MadTest.in 
	$(SETUPTHEPEG) --exitonerror -L .libs -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo MadTest.in 

MadTest.out: MadTest.run
	time $(RUNTHEPEG) -L .libs -d 1 MadTest.run

%.run: %.in done-all-links
	$(SETUPTHEPEG) --exitonerror -L .libs -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo   $<


%.out: %.run
	time $(RUNTHEPEG) -d 1 $<

check-local:
	$(SETUPTHEPEG) --exitonerror -L .libs -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo TestHerwig.in 
	time $(RUNTHEPEG) -L .libs -d 1 TestHerwig.run
	$(SETUPTHEPEG) --exitonerror -L .libs -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo MadTest.in 
	time $(RUNTHEPEG) -L .libs -d 1 MadTest.run

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../TestHerwig.in; $(DESTDIR)$(bindir)/runThePEG TestHerwig.run
	rm -rf _test
	$(mkdir_p) _test
	cd _test; $(LN_S) ../pp2mumu.dat.gz .; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../MadTest.in; $(DESTDIR)$(bindir)/runThePEG MadTest.run
	rm -rf _test


