## workaround for OS X, automake doesn't realize it's meant to use c++
LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
	$(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

LINKDIRS = 

EXTRA_DIST = 

pkglib_LTLIBRARIES = libTheP8I.la

CLEANFILES = 

YODA_LIB = /usr/local/new/lib/libYODA.la

# Version info should be updated if any interface or persistent I/O
# function is changed
libTheP8I_la_LDFLAGS = -version-info 0:0:0 -export-dynamic $(PYTHIA8_LDFLAGS)
libTheP8I_la_SOURCES =
libTheP8I_la_LIBADD = \
	$(top_builddir)/Config/libTheP8IConfig.la \
	$(top_builddir)/Strategy/libTheP8IStrategy.la \
	$(top_builddir)/Hadronization/libTheP8IString.la \
	$(top_builddir)/BE/libTheP8IBE.la $(PYTHIA8_LIBS)

all-local: done-all-links

done-all-links: libTheP8I.la
@EMPTY@ifdef SHOWCOMMAND
	for dir in $(LINKDIRS); do \
	  for sofile in ../$$dir/.libs/*.so* ../$$dir/.libs/*.la; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done; done
	for sofile in .libs/*.so*; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done
	echo "stamp" > done-all-links
@EMPTY@else
	@echo "Sym-linking all shared libraries ..."
	@for dir in $(LINKDIRS); do \
	  for sofile in ../$$dir/.libs/*.so* ../$$dir/.libs/*.la; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done; done
	@for sofile in .libs/*.so*; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done
	@echo "stamp" > done-all-links
@EMPTY@endif

clean-local:
	rm -f *.so* *.la *.rpo done-all-links
