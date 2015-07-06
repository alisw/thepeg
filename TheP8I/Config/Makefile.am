mySOURCES = RndmEngine.cc Pythia8Interface.cc

DOCFILES = TheP8I.h RndmEngine.h Pythia8Interface.h RopeUserHooks.h

INCLUDEFILES = $(DOCFILES)

noinst_LTLIBRARIES = libTheP8IConfig.la
libTheP8IConfig_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

Pythia8Interface.cc: interfaces.pl

include $(top_srcdir)/Config/Makefile.aminclude

