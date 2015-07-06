mySOURCES = ThePEG.cc Containers.cc

DOCFILES = ThePEG.h Constants.h Containers.h Pointers.h \
	   Unitsystem.h algorithm.h \
           std.h Complex.h TemplateTools.h \
	   PhysicalQty.h PhysicalQtyOps.h PhysicalQtyComplex.h

if HAVE_HEPMC
DOCFILES += HepMCHelper.h
endif

INCLUDEFILES = $(DOCFILES)

noinst_LTLIBRARIES = libThePEGConfig.la
libThePEGConfig_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

