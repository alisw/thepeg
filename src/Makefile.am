AUTOMAKE_OPTIONS = -Wno-portability

bin_PROGRAMS = setupThePEG runThePEG
EXTRA_PROGRAMS = runEventLoop

EXTRA_DIST = testpdfs

myLDADD = $(top_builddir)/lib/libThePEG.la 
myLDFLAGS = -export-dynamic

runThePEG_SOURCES = runThePEG.cc
runThePEG_LDADD = $(myLDADD) $(GSLLIBS)
runThePEG_LDFLAGS = $(AM_LDFLAGS) $(myLDFLAGS)

runEventLoop_SOURCES = runEventLoop.cc
runEventLoop_LDADD = -lHepMC $(myLDADD) $(GSLLIBS)
runEventLoop_LDFLAGS = $(AM_LDFLAGS) $(myLDFLAGS)

setupThePEG_SOURCES = setupThePEG.cc
setupThePEG_LDADD = $(myLDADD) $(GSLLIBS)
setupThePEG_LDFLAGS = $(AM_LDFLAGS) $(myLDFLAGS)
setupThePEG_CPPFLAGS = $(AM_CPPFLAGS) \
-DTHEPEG_PKGLIBDIR="\"$(pkglibdir)\""

if USELHAPDF
pkglib_LTLIBRARIES = TestLHAPDF.la
TestLHAPDF_la_SOURCES = TestLHAPDF.cc TestLHAPDF.h
TestLHAPDF_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
endif

dist_pkgdata_DATA = SimpleLEP.in ThePEGDefaults.in ThePEGParticles.in debugItems.txt TestLHAPDF.in MultiLEP.in

rpodir = $(pkglibdir)
nodist_rpo_DATA = ThePEGDefaults.rpo

CLEANFILES = SimpleLEP.log SimpleLEP.out SimpleLEP.run SimpleLEP.tex \
             MultiLEP.log MultiLEP.out MultiLEP.run MultiLEP.tex \
             ThePEGDefaults.rpo .done-all-links \
             TestLHAPDF.log TestLHAPDF.out TestLHAPDF.run TestLHAPDF.tex \
             .runThePEG.timer.TestLHAPDF.run SimpleLEP.dump MultiLEP.dump 

save:
	mkdir -p save
	cp *.log *.out *.run *.tex *.rpo save

valgrind:
	valgrind --leak-check=full --num-callers=25 --track-fds=yes --freelist-vol=100000000 --leak-resolution=med --trace-children=yes ./setupThePEG -r ThePEGDefaults.rpo SimpleLEP.in &> /tmp/valgrind.out
	valgrind --leak-check=full --num-callers=25 --track-fds=yes --freelist-vol=100000000 --leak-resolution=med --trace-children=yes ./runThePEG SimpleLEP.run >> /tmp/valgrind.out 2>&1

INPUTFILES = ThePEGDefaults.in ThePEGParticles.in \
             SimpleLEP.in SimpleLEP.mod MultiLEP.in TestLHAPDF.in

.done-all-links:
@EMPTY@ifdef SHOWCOMMAND
	for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	echo "stamp" > .done-all-links
@EMPTY@else
	@echo "sym-linking input files files..."
	@for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	@echo "stamp" > .done-all-links
@EMPTY@endif

ThePEGDefaults.rpo: .done-all-links setupThePEG ThePEGDefaults.in ThePEGParticles.in
	./setupThePEG --init --exitonerror -L../lib -r ThePEGDefaults.rpo ThePEGDefaults.in
	cd ../lib; rm -f ThePEGDefaults.rpo; $(LN_S) -f ../src/ThePEGDefaults.rpo .

check-local: ThePEGDefaults.rpo
	./setupThePEG --exitonerror -r ThePEGDefaults.rpo SimpleLEP.in
	mv SimpleLEP.out SimpleLEP.cmp
	time ./runThePEG -d 0 SimpleLEP.run
	diff SimpleLEP.out SimpleLEP.cmp
	mv SimpleLEP.out SimpleLEP.cmp
	time ./runThePEG --resume -d 0 SimpleLEP.dump
	diff SimpleLEP.out SimpleLEP.cmp
	rm SimpleLEP.cmp
	time ./runThePEG -d 0 -m SimpleLEP.mod SimpleLEP.run
	./setupThePEG --exitonerror -r ThePEGDefaults.rpo MultiLEP.in
	time ./runThePEG -d 0 MultiLEP.run
if USELHAPDF
	LHAPATH=$(srcdir)/testpdfs ./setupThePEG --exitonerror -r ThePEGDefaults.rpo TestLHAPDF.in
	LHAPATH=$(srcdir)/testpdfs time ./runThePEG -d 1 -x .libs/TestLHAPDF.so TestLHAPDF.run
endif

SimpleLEP.run: .done-all-links setupThePEG ThePEGDefaults.rpo SimpleLEP.in
	./setupThePEG --exitonerror -r ThePEGDefaults.rpo SimpleLEP.in

SimpleLEP.out: runThePEG SimpleLEP.run
	time ./runThePEG -d 0 SimpleLEP.run

MultiLEP.run: .done-all-links setupThePEG ThePEGDefaults.rpo MultiLEP.in
	./setupThePEG --exitonerror -r ThePEGDefaults.rpo MultiLEP.in

MultiLEP.out: runThePEG MultiLEP.run
	time ./runThePEG -d 0 MultiLEP.run

if USELHAPDF
TestLHAPDF.run: .done-all-links setupThePEG ThePEGDefaults.rpo TestLHAPDF.in TestLHAPDF.la
	LHAPATH=$(srcdir)/testpdfs ./setupThePEG --exitonerror -r ThePEGDefaults.rpo TestLHAPDF.in

TestLHAPDF.out: runThePEG TestLHAPDF.run
	LHAPATH=$(srcdir)/testpdfs time ./runThePEG -d 1 -x .libs/TestLHAPDF.so TestLHAPDF.run
endif

install-exec-hook:
	cd $(DESTDIR)$(pkglibdir); \
	rm -f setupThePEG-$(VERSION) runThePEG-$(VERSION); \
	cp $(DESTDIR)$(bindir)/setupThePEG setupThePEG-$(VERSION); \
	cp $(DESTDIR)$(bindir)/runThePEG runThePEG-$(VERSION)

install-data-hook:
	rm -f $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo
	./setupThePEG --init --exitonerror -L$(DESTDIR)$(pkglibdir) -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(DESTDIR)$(pkgdatadir)/ThePEGDefaults.in
if USELHAPDF
#	LHAPATH=$(srcdir)/testpdfs ./setupThePEG --init --exitonerror -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo
endif
	cd $(DESTDIR)$(pkglibdir); \
	if test -n "$(DESTDIR)"; \
	then sed -i.bak -e "s@$(DESTDIR)@@g" ThePEGDefaults.rpo; \
	rm -f ThePEGDefaults.rpo.bak; \
	fi ; \
	cp ThePEGDefaults.rpo ThePEGDefaults-$(VERSION).rpo

uninstall-local:
	cd $(DESTDIR)$(pkglibdir); \
	rm -f setupThePEG-$(VERSION) runThePEG-$(VERSION) ThePEGDefaults-$(VERSION).rpo 

installcheck-local:
	$(mkdir_p) _test_thepeg_installation_
	cd _test_thepeg_installation_; \
$(DESTDIR)$(bindir)/setupThePEG --exitonerror $(DESTDIR)$(pkgdatadir)/SimpleLEP.in; \
mv SimpleLEP.out SimpleLEP.cmp; \
$(DESTDIR)$(bindir)/runThePEG SimpleLEP.run; \
diff SimpleLEP.out SimpleLEP.cmp; \
mv SimpleLEP.out SimpleLEP.cmp; \
$(DESTDIR)$(bindir)/runThePEG --resume SimpleLEP.dump; \
diff SimpleLEP.out SimpleLEP.cmp; \
$(DESTDIR)$(bindir)/setupThePEG --exitonerror $(DESTDIR)$(pkgdatadir)/MultiLEP.in; \
$(DESTDIR)$(bindir)/runThePEG MultiLEP.run;
if USELHAPDF
	export LHAPATH=$$(pwd)/$(srcdir)/testpdfs; cd _test_thepeg_installation_; \
$(DESTDIR)$(bindir)/setupThePEG --exitonerror $(DESTDIR)$(pkgdatadir)/TestLHAPDF.in; \
$(DESTDIR)$(bindir)/runThePEG -x $(DESTDIR)$(pkglibdir)/TestLHAPDF.so TestLHAPDF.run
endif
	rm -rf _test_thepeg_installation_

%.run: %.in .done-all-links ThePEGDefaults.rpo
	./setupThePEG --exitonerror -L .libs -r ../lib/ThePEGDefaults.rpo $<

%.out: %.run
	time ./runThePEG -d 0 $<

