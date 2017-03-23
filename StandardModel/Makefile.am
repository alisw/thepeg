mySOURCES = StandardModelBase.cc RunningCoupling.cc AlphaEMBase.cc CKMBase.cc AlphaSBase.cc

DOCFILES = RunningCoupling.h AlphaEMBase.h AlphaSBase.h CKMBase.h  \
           StandardModelBase.h

INCLUDEFILES = $(DOCFILES) StandardModelBase.fh


noinst_LTLIBRARIES = libThePEGSM.la
pkglib_LTLIBRARIES = O1AlphaS.la SimpleAlphaEM.la StandardCKM.la

libThePEGSM_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

# Version info should be updated if any interface or persistent I/O
# function is changed
O1AlphaS_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
O1AlphaS_la_SOURCES = O1AlphaS.cc O1AlphaS.h

# Version info should be updated if any interface or persistent I/O
# function is changed
SimpleAlphaEM_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
SimpleAlphaEM_la_SOURCES = SimpleAlphaEM.cc SimpleAlphaEM.h

# Version info should be updated if any interface or persistent I/O
# function is changed
StandardCKM_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
StandardCKM_la_SOURCES = StandardCKM.cc StandardCKM.h

include $(top_srcdir)/Config/Makefile.aminclude

