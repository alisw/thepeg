mySOURCES = Amplitude.cc ColourLines.cc DiagramBase.cc Tree2toNDiagram.cc BlobDiagram.cc \
            MEBase.cc BlobMEBase.cc MEGroup.cc ReweightBase.cc ME2to2Base.cc ME2to2QCD.cc

myQCDSOURCES = 

DOCFILES = Amplitude.h ColourLines.h DiagramBase.h Tree2toNDiagram.h BlobDiagram.h MEBase.h \
           BlobMEBase.h MEGroup.h ReweightBase.h ME2to2Base.h ME2to2QCD.h

INCLUDEFILES = $(DOCFILES) ColourLines.xh \
               DiagramBase.xh DiagramBase.fh \
               Tree2toNDiagram.xh MEBase.fh MEGroup.fh

noinst_LTLIBRARIES = libThePEGME.la
pkglib_LTLIBRARIES = MEQCD.la MEee2gZ2qq.la MENCDIS.la \
                     ReweightConstant.la ReweightMinPT.la

libThePEGME_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

# Version info should be updated if any interface or persistent I/O
# function is changed
MEQCD_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
MEQCD_la_SOURCES = MEGG2GG.cc MEQG2QG.cc MEQQ2QQ3.cc MEQq2Qq4.cc MEGG2QQ.cc \
                   MEQQ2GG.cc MEQQ2qq2.cc MEqq2qq1.cc MEGG2GG.h MEQG2QG.h \
                   MEQQ2QQ3.h MEQq2Qq4.h MEGG2QQ.h MEQQ2GG.h MEQQ2qq2.h \
                   MEqq2qq1.h

# Version info should be updated if any interface or persistent I/O
# function is changed
MEee2gZ2qq_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
MEee2gZ2qq_la_SOURCES = MEee2gZ2qq.cc MEee2gZ2qq.h

# Version info should be updated if any interface or persistent I/O
# function is changed
ReweightConstant_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
ReweightConstant_la_SOURCES = ReweightConstant.cc ReweightConstant.h

# Version info should be updated if any interface or persistent I/O
# function is changed
ReweightMinPT_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
ReweightMinPT_la_SOURCES = ReweightMinPT.cc ReweightMinPT.h

# Version info should be updated if any interface or persistent I/O
# function is changed
MENCDIS_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
MENCDIS_la_SOURCES = MENCDIS.cc MENCDIS.h

include $(top_srcdir)/Config/Makefile.aminclude

