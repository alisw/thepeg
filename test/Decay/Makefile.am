CLEANFILES = TestDecays.log TestDecays.out TestDecays.run TestDecays.tex \
             done-all-links .runThePEG.timer.TestDecays.run

pkglib_LTLIBRARIES = TestDecays.la

# Version info should be updated if any interface or persistent I/O
# function is changed
TestDecays_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
TestDecays_la_SOURCES = TestDecays.cc TestDecays.h TestDecays.icc

EXTRA_DIST = TestDecays.in

INPUTFILES = TestDecays.in

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

TestDecays.run: done-all-links TestDecays.in
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo TestDecays.in

TestDecays.out: TestDecays.run TestDecays.la
	time $(RUNTHEPEG) -d 1 -x .libs/TestDecays.so TestDecays.run

check-local:
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo TestDecays.in
	time $(RUNTHEPEG) -d 1 -x .libs/TestDecays.so TestDecays.run

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../TestDecays.in; $(DESTDIR)$(bindir)/runThePEG -x $(DESTDIR)$(pkglibdir)/TestDecays.so TestDecays.run
	rm -rf _test

