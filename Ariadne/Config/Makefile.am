mySOURCES = CloneBase.cc Settings.cc

DOCFILES = Ariadne5.h CloneBase.h UnitFO.h Settings.h

INCLUDEFILES = $(DOCFILES) config.h CloneBase.fh CloneBase.icc

EXTRA_DIST = $(INCLUDEFILES)

noinst_LTLIBRARIES = libAriadne5Config.la

libAriadne5Config_la_SOURCES =  $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude
