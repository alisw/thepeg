AUTOMAKE_OPTIONS = -Wno-portability

dist_pkgdata_DATA = Ariadne5Defaults.in Ariadne5Remove.in TestAriadne5.in \
                    TestAriadne5DIS.in TestAriadne5Z0.in 

rpodir = $(pkglibdir)
nodist_rpo_DATA = Ariadne5Defaults.rpo

CLEANFILES = Ariadne5Defaults.rpo done-all-links \
             TestAriadne5.log TestAriadne5.out TestAriadne5.run TestAriadne5.tex \
             TestAriadne5Z0.log TestAriadne5Z0.out TestAriadne5Z0.run \
	     TestAriadne5Z0.tex TestAriadne5DIS.log TestAriadne5DIS.out \
             TestAriadne5DIS.run TestAriadne5DIS.tex

INPUTFILES = Ariadne5Defaults.in Ariadne5Remove.in TestAriadne5.in \
             TestAriadne5Z0.in TestAriadne5DIS.in

pkglib_LTLIBRARIES = MEqq2gZ2ll.la DYAnalysis.la DISAnalysis.la ETFlow.la

# Version info should be updated if any interface or persistent I/O
# function is changed
MEqq2gZ2ll_la_LDFLAGS = -module -version-info 1:0:0
MEqq2gZ2ll_la_SOURCES = MEqq2gZ2ll.cc MEqq2gZ2ll.h MEqq2gZ2ll.icc

# Version info should be updated if any interface or persistent I/O
# function is changed
DYAnalysis_la_LDFLAGS = -module -version-info 1:0:0
DYAnalysis_la_SOURCES = DYAnalysis.cc DYAnalysis.h DYAnalysis.icc

DISAnalysis_la_LDFLAGS = -module -version-info 1:0:0
DISAnalysis_la_SOURCES = DISAnalysis.cc DISAnalysis.h DISAnalysis.icc

ETFlow_la_LDFLAGS = -module -version-info 1:0:0
ETFlow_la_SOURCES = ETFlow.cc ETFlow.h

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

Ariadne5Defaults.rpo: done-all-links Ariadne5Defaults.in ../lib/libAriadne5.so $(THEPEGLIB)/ThePEGDefaults.rpo  $(THEP8ILIBPATH)/TheP8IDefaults.rpo MEqq2gZ2ll.la
	$(SETUPTHEPEG) --exitonerror --init -L../lib -L$(THEP8ILIBPATH) -L$(THEPEGLIB) -L.libs -r $(THEP8ILIBPATH)/TheP8IDefaults.rpo -o Ariadne5Defaults.rpo Ariadne5Defaults.in
	cd ../lib; $(LN_S) -f ../src/Ariadne5Defaults.rpo .

check-local:
	$(SETUPTHEPEG) --exitonerror -r ../lib/Ariadne5Defaults.rpo TestAriadne5.in
	time $(RUNTHEPEG) -d 0 TestAriadne5.run
	$(SETUPTHEPEG) --exitonerror -r ../lib/Ariadne5Defaults.rpo TestAriadne5DIS.in
	time $(RUNTHEPEG) -d 0 TestAriadne5DIS.run
	$(SETUPTHEPEG) --exitonerror -r ../lib/Ariadne5Defaults.rpo TestAriadne5Z0.in
	time $(RUNTHEPEG) -d 0 TestAriadne5Z0.run

%.run: %.in done-all-links Ariadne5Defaults.rpo
	$(SETUPTHEPEG) --exitonerror -r ../lib/Ariadne5Defaults.rpo $<

%.out: %.run
	time $(RUNTHEPEG) $<

install-data-local:
	LD_LIBRARY_PATH=$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH $(DESTDIR)$(bindir)/setupThePEG --exitonerror --init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(srcdir)/Ariadne5Defaults.in

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(DESTDIR)$(bindir)/setupThePEG $(DESTDIR)$(pkgdatadir)/TestAriadne5.in; $(DESTDIR)$(bindir)/runThePEG TestAriadne5.run; $(DESTDIR)$(bindir)/setupThePEG $(DESTDIR)$(pkgdatadir)/TestAriadne5DIS.in; $(DESTDIR)$(bindir)/runThePEG TestAriadne5DIS.run; $(DESTDIR)$(bindir)/setupThePEG $(DESTDIR)$(pkgdatadir)/TestAriadne5Z0.in; $(DESTDIR)$(bindir)/runThePEG TestAriadne5Z0.run
	rm -rf _test

save:
	mkdir -p save
	cp *.rpo *.out *.log *.tex .*.timer.* save
