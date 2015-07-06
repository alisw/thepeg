mySOURCES = String.cc StringRegion.cc EndPoint.cc Oriented.cc \
	    LundZGenerator.cc LundPtGenerator.cc \
	    LundFlavourGenerator.cc LundFragHandler.cc LundFlavourHandler.cc \
            StringCollapser.cc

DOCFILES = LundFragHandler.h FragConfig.h Oriented.h \
           OrientedIndex.h EndPoint.h \
           StringRegion.h String.h Hadron.h LundFlavourHandler.h \
           LundZGenerator.h LundPtGenerator.h LundFlavourGenerator.h \
           StringCollapser.h

INCLUDEFILES = $(DOCFILES)  LundFragHandler.xh Oriented.icc \
               String.icc StringRegion.icc EndPoint.icc \
               Hadron.icc LundFragHandler.icc  LundPtGenerator.icc \
               LundZGenerator.icc LundFlavourGenerator.icc \
               LundFlavourHandler.icc StringCollapser.icc

pkglib_LTLIBRARIES = libP7String.la

# Version info should be updated if any interface or persistent I/O
# function is changed
libP7String_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
libP7String_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

