CLEANFILES = TestMisc.log TestMisc.out TestMisc.run TestMisc.tex \
             done-all-links .runThePEG.timer.TestMisc.run

pkglib_LTLIBRARIES = TestMisc.la

# Version info should be updated if any interface or persistent I/O
# function is changed
TestMisc_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
TestMisc_la_SOURCES = TestMisc.cc TestMisc.h TestMisc.icc

EXTRA_DIST = TestMisc.in

INPUTFILES = TestMisc.in

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

TestMisc.run: done-all-links TestMisc.in
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo TestMisc.in

TestMisc.out: TestMisc.run TestMisc.la
	time $(RUNTHEPEG) -d 1 -x .libs/TestMisc.so TestMisc.run

check-local:
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo TestMisc.in
	time $(RUNTHEPEG) -d 1 -x .libs/TestMisc.so TestMisc.run

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../TestMisc.in; $(DESTDIR)$(bindir)/runThePEG -x $(DESTDIR)$(pkglibdir)/TestMisc.so TestMisc.run
	rm -rf _test

