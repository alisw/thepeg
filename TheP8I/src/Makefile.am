
dist_pkgdata_DATA = SimpleLEP.in SimpleLHC.in TheP8IDefaults.in \
                    TheP8IRemove.in TheP8IParticles.in TestJunctions.in \
                    OverlapStringFragmenter.in

CLEANFILES = SimpleLEP.log SimpleLEP.out SimpleLEP.run SimpleLEP.tex SimpleLEP.dump \
             SimpleLHC.log SimpleLHC.out SimpleLHC.run SimpleLHC.tex \
             TestJunctions.log TestJunctions.out TestJunctions.run TestJunctions.tex \
             TheP8IDefaults.rpo TheP8IParticles.in done-all-links

pkglib_LTLIBRARIES = TestJunctions.la

rpodir = $(pkglibdir)
nodist_rpo_DATA = TheP8IDefaults.rpo

# Version info should be updated if any interface or persistent I/O
# function is changed
TestJunctions_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
TestJunctions_la_SOURCES = TestJunctions.cc TestJunctions.h TestJunctions.icc

INPUTFILES = TheP8IDefaults.in TheP8IRemove.in SimpleLEP.in SimpleLHC.in \
             TestJunctions.in OverlapStringFragmenter.in

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

TheP8IParticles.in: $(builddir)/P8Particles.pl
	$(builddir)/P8Particles.pl @PYTHIA8_DIR@/xmldoc/ParticleData.xml > TheP8IParticles.in

TheP8IDefaults.rpo: done-all-links TheP8IDefaults.in TheP8IParticles.in ../lib/libTheP8I.so  $(THEPEGLIB)/ThePEGDefaults.rpo
	$(SETUPTHEPEG) -L../lib -L$(THEPEGLIB) --exitonerror --init -r $(THEPEGLIB)/ThePEGDefaults.rpo -o TheP8IDefaults.rpo TheP8IDefaults.in
	cd ../lib; rm -f TheP8IDefaults.rpo; $(LN_S) -f ../src/TheP8IDefaults.rpo .

check-local: TheP8IDefaults.rpo
	$(SETUPTHEPEG) --exitonerror -r ../lib/TheP8IDefaults.rpo SimpleLEP.in
	time $(RUNTHEPEG) -d 0 SimpleLEP.run
	$(SETUPTHEPEG) --exitonerror -r ../lib/TheP8IDefaults.rpo TestJunctions.in
	time $(RUNTHEPEG) -d 1 -x .libs/TestJunctions.so TestJunctions.run
#	$(SETUPTHEPEG) --exitonerror -r ../lib/TheP8IDefaults.rpo SimpleLHC.in
#	time $(RUNTHEPEG) -d 1 SimpleLHC.run

SimpleLEP.run: done-all-links SimpleLEP.in TheP8IDefaults.rpo
	$(SETUPTHEPEG) --exitonerror -r ../lib/TheP8IDefaults.rpo SimpleLEP.in

SimpleLHC.run: done-all-links SimpleLHC.in TheP8IDefaults.rpo
	$(SETUPTHEPEG) --exitonerror -r ../lib/TheP8IDefaults.rpo SimpleLHC.in

SimpleLEP.out: SimpleLEP.run
	time $(RUNTHEPEG) -d 0 SimpleLEP.run

SimpleLHC.out: SimpleLHC.run
	time $(RUNTHEPEG) -d 1 SimpleLHC.run

TestJunctions.run: done-all-links TestJunctions.in
	$(SETUPTHEPEG) --exitonerror -r ../lib/TheP8IDefaults.rpo TestJunctions.in

TestJunctions.out: TestJunctions.run TestJunctions.la
	time $(RUNTHEPEG) -d 1 -x .libs/TestJunctions.so TestJunctions.run

install-data-local:
	LD_LIBRARY_PATH=$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH $(DESTDIR)$(bindir)/setupThePEG --exitonerror --init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(srcdir)/TheP8IDefaults.in

unregister-from-repo:
	LD_LIBRARY_PATH=$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH $(DESTDIR)$(bindir)/setupThePEG --exitonerror --init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(srcdir)/TheP8IRemove.in

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(DESTDIR)$(bindir)/setupThePEG --exitonerror $(DESTDIR)$(pkgdatadir)/SimpleLEP.in; $(DESTDIR)$(bindir)/runThePEG SimpleLEP.run; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../TestJunctions.in; $(DESTDIR)$(bindir)/runThePEG -x $(DESTDIR)$(pkglibdir)/TestJunctions.so TestJunctions.run
	rm -rf _test

save:
	mkdir -p save
	cp *.log *.out *.run *.tex *.rpo save

