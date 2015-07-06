mySOURCES = Pythia7Strategy.cc

DOCFILES = Pythia7Strategy.h

INCLUDEFILES = $(DOCFILES) Pythia7Strategy.icc

pkglib_LTLIBRARIES = Pythia7Strategy.la

# Version info should be updated if any interface or persistent I/O
# function is changed
Pythia7Strategy_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
Pythia7Strategy_la_SOURCES = Pythia7Strategy.cc Pythia7Strategy.h \
                             Pythia7Strategy.icc

include $(top_srcdir)/Config/Makefile.aminclude

