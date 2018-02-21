mySOURCES = EventGenerator.cc RandomGenerator.cc Strategy.cc \
            BaseRepository.cc Repository.cc StandardRandom.cc  \
            UseRandom.cc CurrentGenerator.cc Main.cc

DOCFILES = BaseRepository.h EventGenerator.h RandomGenerator.h \
           Repository.h StandardRandom.h Strategy.h  \
           UseRandom.h CurrentGenerator.h Main.h

INCLUDEFILES = $(DOCFILES) BaseRepository.tcc \
               BaseRepository.xh EventGenerator.fh \
               EventGenerator.tcc EventGenerator.xh \
               Strategy.fh RandomGenerator.fh UseRandom.fh \
               CurrentGenerator.fh

noinst_LTLIBRARIES = libThePEGRepo.la
pkglib_LTLIBRARIES = ThePEGStrategy.la MultiEventGenerator.la

libThePEGRepo_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
libThePEGRepo_la_CPPFLAGS = $(AM_CPPFLAGS) \
-DTHEPEG_PKGLIBDIR="\"$(pkglibdir)\"" \
-DTHEPEG_PKGDATADIR="\"$(pkgdatadir)\""

nodist_libThePEGRepo_la_SOURCES = versionstamp.inc
BUILT_SOURCES = versionstamp.inc
CLEANFILES = versionstamp.inc

AUTOMAKE_OPTIONS = -Wno-portability

VERSIONSTAMP := $(shell hg -R $(top_srcdir) parents --template '"ThePEG {node|short} ({branch})"' 2> /dev/null || echo \"$(PACKAGE_STRING)\" || true )

.PHONY: update_versionstamp
versionstamp.inc: update_versionstamp
	@[ -f $@ ] || touch $@
	@echo '$(VERSIONSTAMP)' | cmp -s $@ - || echo '$(VERSIONSTAMP)' > $@


# Version info should be updated if any interface or persistent I/O
# function is changed
ThePEGStrategy_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
ThePEGStrategy_la_SOURCES = ThePEGStrategy.cc  ThePEGStrategy.h

# Version info should be updated if any interface or persistent I/O
# function is changed
MultiEventGenerator_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
MultiEventGenerator_la_SOURCES = MultiEventGenerator.cc MultiEventGenerator.h

include $(top_srcdir)/Config/Makefile.aminclude

# Compile and use Boost unit tests only if boost unit test libs are available
check_PROGRAMS =
repository_test_SOURCES =
repository_test_LDADD =
repository_test_LDFLAGS =
repository_test_CPPFLAGS =
TESTS =

if COND_BOOSTTEST
 THEPEGLDADD = $(top_builddir)/lib/libThePEG.la 
 check_PROGRAMS += repository_test
 repository_test_SOURCES += tests/repositoryTestsMain.cc \
 tests/repositoryTestsGlobalFixture.h \
 tests/repositoryTestRandomGenerator.h
 repository_test_LDADD += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) $(THEPEGLDADD) $(GSLLIBS) 
 repository_test_LDFLAGS += $(AM_LDFLAGS) -export-dynamic $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) 
 repository_test_CPPFLAGS += $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -DTHEPEG_PKGLIBDIR="\"$(pkglibdir)\"" -DTHEPEG_PKGDATADIR="\"$(pkgdatadir)\""
 TESTS += repository_test
endif
