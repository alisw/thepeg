
dist_pkgdata_DATA = SimpleLEP.in RivetLEP.in SimpleLHC.in Pythia7Defaults.in \
                    Pythia7Remove.in

rpodir = $(pkglibdir)
nodist_rpo_DATA = Pythia7Defaults.rpo

CLEANFILES = SimpleLEP.log SimpleLEP.out SimpleLEP.run SimpleLEP.tex \
             RivetLEP.log RivetLEP.out RivetLEP.run RivetLEP.tex \
             SimpleLHC.log SimpleLHC.out SimpleLHC.run SimpleLHC.tex \
             RivetLEP.aida RivetLEP.yoda Pythia7Defaults.rpo SimpleLHC.aida \
             SimpleLHC.dat done-all-links

INPUTFILES = Pythia7Defaults.in Pythia7Remove.in SimpleLEP.in RivetLEP.in \
             SimpleLHC.in

pkglib_LTLIBRARIES = CMSlrnsac.la
CMSlrnsac_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
CMSlrnsac_la_SOURCES = CMSlrnsac.cc CMSlrnsac.h


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

Pythia7Defaults.rpo: done-all-links Pythia7Defaults.in ../lib/Pythia7Strategy.so ../lib/libP7String.so ../lib/libP7Shower.so $(THEPEGLIB)/ThePEGDefaults.rpo
	$(SETUPTHEPEG) -L../lib -L$(THEPEGLIB) --exitonerror --init -r $(THEPEGLIB)/ThePEGDefaults.rpo -o Pythia7Defaults.rpo Pythia7Defaults.in
	cd ../lib; rm -f Pythia7Defaults.rpo; $(LN_S) -f ../src/Pythia7Defaults.rpo .

check-local: Pythia7Defaults.rpo CMSlrnsac.la
	$(SETUPTHEPEG) --exitonerror -r ../lib/Pythia7Defaults.rpo SimpleLEP.in
	time $(RUNTHEPEG) -d 0 SimpleLEP.run
	$(SETUPTHEPEG) --exitonerror -L.libs -r ../lib/Pythia7Defaults.rpo SimpleLHC.in
	time $(RUNTHEPEG) -d 1 SimpleLHC.run
if HAVE_RIVET
	$(SETUPTHEPEG) --exitonerror -r ../lib/Pythia7Defaults.rpo RivetLEP.in
	time $(RUNTHEPEG) -d 0 RivetLEP.run
endif

SimpleLEP.run: done-all-links SimpleLEP.in Pythia7Defaults.rpo
	$(SETUPTHEPEG) -L.lib --exitonerror -r ../lib/Pythia7Defaults.rpo SimpleLEP.in

SimpleLEP.out: SimpleLEP.run
	time $(RUNTHEPEG) -d 0 SimpleLEP.run

SimpleLHC.run: CMSlrnsac.la done-all-links SimpleLHC.in Pythia7Defaults.rpo
	$(SETUPTHEPEG) --exitonerror -L.libs -r ../lib/Pythia7Defaults.rpo SimpleLHC.in

SimpleLHC.out: SimpleLHC.run
	time $(RUNTHEPEG) -d 1 SimpleLHC.run

if HAVE_RIVET
RivetLEP.run: done-all-links RivetLEP.in Pythia7Defaults.rpo
	$(SETUPTHEPEG) --exitonerror -r ../lib/Pythia7Defaults.rpo RivetLEP.in

RivetLEP.out: RivetLEP.run
	time $(RUNTHEPEG) -d 0 RivetLEP.run
endif

install-data-local:
	LD_LIBRARY_PATH=$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH $(DESTDIR)$(bindir)/setupThePEG --exitonerror --init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(srcdir)/Pythia7Defaults.in

unregister-from-repo:
	LD_LIBRARY_PATH=$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH $(DESTDIR)$(bindir)/setupThePEG --exitonerror --init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(srcdir)/Pythia7Remove.in

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(DESTDIR)$(bindir)/setupThePEG --exitonerror $(DESTDIR)$(pkgdatadir)/SimpleLEP.in; $(DESTDIR)$(bindir)/runThePEG SimpleLEP.run; $(DESTDIR)$(bindir)/setupThePEG --exitonerror $(DESTDIR)$(pkgdatadir)/SimpleLHC.in; $(DESTDIR)$(bindir)/runThePEG SimpleLHC.run
if HAVE_RIVET
	cd _test; $(DESTDIR)$(bindir)/setupThePEG --exitonerror $(DESTDIR)$(pkgdatadir)/RivetLEP.in; $(DESTDIR)$(bindir)/runThePEG RivetLEP.run
endif
	rm -rf _test

save:
	mkdir -p save
	cp *.log *.out *.run *.tex *.rpo save

