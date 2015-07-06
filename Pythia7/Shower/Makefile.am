mySOURCES = Basics.cc Beam.cc Shower.cc TimeShower.cc SpaceShower.cc \
            ShowerHandler.cc SpaceShowerHandler.cc TimeShowerHandler.cc

DOCFILES = Basics.h Beam.h Shower.h TimeShower.h SpaceShower.h \
           ShowerHandler.h SpaceShowerHandler.h TimeShowerHandler.h

INCLUDEFILES = $(DOCFILES) ShowerHandler.icc SpaceShowerHandler.icc \
               TimeShowerHandler.icc SpaceShowerHandler.fh TimeShowerHandler.fh

pkglib_LTLIBRARIES = libP7Shower.la

# Version info should be updated if any interface or persistent I/O
# function is changed
libP7Shower_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
libP7Shower_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

