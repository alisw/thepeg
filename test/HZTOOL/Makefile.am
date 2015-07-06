CLEANFILES = TestHZTool.log TestHZTool.out TestHZTool.run TestHZTool.tex \
             testhztool.rz TestHZTool.hzlog done-all-links \
             .runThePEG.timer.TestHZTool.run TestHZTool-Plotter-1.dot

if USEHZTOOL
if USEHEPMC

EXTRA_DIST = TestHZTool.in

INPUTFILES = TestHZTool.in

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

pkglib_LTLIBRARIES = HZTool.la

# Version info should be updated if any interface or persistent I/O
# function is changed
HZTool_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO) $(HZTOOL_LDFLAGS) $(HEPMC_LDFLAGS)
HZTool_la_SOURCES = HZTool.cc HZTool.h HZTool.icc hztoolinit.f hz98051.f
HZTool_la_LIBADD = $(HZTOOL_LIBS) $(HEPMC_LIBS)

check-local:
	$(SETUPTHEPEG) --exitonerror -r $(PYTHIA7LIBPATH)/Pythia7Defaults.rpo -L .libs TestHZTool.in
	time $(RUNTHEPEG) TestHZTool.run

installcheck-local:
	$(mkdir_p) _test
	cd _test; $(DESTDIR)$(bindir)/setupThePEG --exitonerror ../TestHZTool.in; $(DESTDIR)$(bindir)/runThePEG TestHZTool.run
	rm -rf _test

endif
endif


