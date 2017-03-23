mySOURCES = PersistentOStream.cc PersistentIStream.cc

DOCFILES = InputDescription.h PersistentIStream.h PersistentOStream.h

INCLUDEFILES = $(DOCFILES) PersistentIStream.fh PersistentIStream.xh \
            PersistentOStream.fh PersistentOStream.xh 

noinst_LTLIBRARIES = libThePEGPersistency.la
libThePEGPersistency_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

