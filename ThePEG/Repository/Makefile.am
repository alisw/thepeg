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


# Version info should be updated if any interface or persistent I/O
# function is changed
ThePEGStrategy_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
ThePEGStrategy_la_SOURCES = ThePEGStrategy.cc  ThePEGStrategy.h


# Version info should be updated if any interface or persistent I/O
# function is changed
MultiEventGenerator_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
MultiEventGenerator_la_SOURCES = MultiEventGenerator.cc MultiEventGenerator.h

include $(top_srcdir)/Config/Makefile.aminclude

