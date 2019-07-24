## workaround for OS X, automake doesn't realize it's meant to use c++
LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
	$(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

LINKDIRS = Config Utilities Helicity Interface LesHouches Vectors \
           PDT PDF Persistency Handlers MatrixElement Pointer ACDC \
           StandardModel Repository EventRecord Cuts Analysis

EXTRA_DIST = Makefile.dist.install

dist_pkgdata_DATA = ThePEG.el

pkglib_LTLIBRARIES = libThePEG.la

CLEANFILES = Makefile.dist Makefile.common

# Version info should be updated if any interface or persistent I/O
# function is changed
libThePEG_la_LDFLAGS = $(AM_LDFLAGS) $(LIBTOOLVERSIONINFO) -export-dynamic
libThePEG_la_SOURCES =
libThePEG_la_LIBADD = \
	$(top_builddir)/Utilities/libThePEGUtilities.la \
	$(top_builddir)/Vectors/libThePEGVectors.la \
	$(top_builddir)/Persistency/libThePEGPersistency.la \
	$(top_builddir)/Pointer/libThePEGReferenceCounted.la \
	$(top_builddir)/Config/libThePEGConfig.la \
	$(top_builddir)/PDF/libThePEGPDF.la \
	$(top_builddir)/MatrixElement/libThePEGME.la \
	$(top_builddir)/Repository/libThePEGRepo.la \
	$(top_builddir)/Handlers/libThePEGHandlers.la \
	$(top_builddir)/EventRecord/libThePEGEvent.la \
	$(top_builddir)/PDT/libThePEGPDT.la \
	$(top_builddir)/StandardModel/libThePEGSM.la \
	$(top_builddir)/Interface/libThePEGInter.la \
	$(top_builddir)/Helicity/libThePEGHelicity.la \
	$(top_builddir)/Cuts/libThePEGCuts.la \
	$(top_builddir)/Analysis/libThePEGHist.la \
	-ldl $(LIBS) \
	$(GSLLIBS)

all-local: .done-all-links

.done-all-links: libThePEG.la
@EMPTY@ifdef SHOWCOMMAND
	for dir in $(LINKDIRS); do \
	  for sofile in ../$$dir/.libs/*.so* ../$$dir/.libs/*.la; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done; done
	for sofile in .libs/*.so*; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done
	echo "stamp" > .done-all-links
@EMPTY@else
	@echo "Sym-linking all shared libraries ..."
	@for dir in $(LINKDIRS); do \
	  for sofile in ../$$dir/.libs/*.so* ../$$dir/.libs/*.la; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done; done
	@for sofile in .libs/*.so*; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done
	@echo "stamp" > .done-all-links
@EMPTY@endif

install-data-local: Makefile.common.install
	sed -e s:@THEPEGLIB@:$(pkglibdir):g -e s:@THEPEGINCLUDE@:-I$(includedir):g Makefile.common.install > Makefile.common
	sed -e s:@THEPEGLIB@:$(pkglibdir):g $(srcdir)/Makefile.dist.install > Makefile.dist
	test -z "$(pkglibdir)" || $(mkdir_p) "$(DESTDIR)$(pkglibdir)"
	$(INSTALL_DATA) Makefile.common $(DESTDIR)$(pkglibdir)/Makefile.common
	$(INSTALL_DATA) Makefile.dist $(DESTDIR)$(pkglibdir)/Makefile

uninstall-local:
	rm -f $(DESTDIR)$(pkglibdir)/Makefile.common $(DESTDIR)$(pkglibdir)/Makefile

clean-local:
	rm -f *.so* *.la *.rpo .done-all-links
