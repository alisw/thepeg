mySOURCES = InterfaceBase.cc InterfacedBase.cc Interfaced.cc ParVector.cc \
          Parameter.cc RefVector.cc Reference.cc Switch.cc Command.cc \
          ClassDocumentation.cc Deleted.cc

DOCFILES = Command.h InterfaceBase.h Interfaced.h InterfacedBase.h ParVector.h \
           Parameter.h RefVector.h Reference.h Switch.h ClassDocumentation.h \
           Interface.h

INCLUDEFILES = $(DOCFILES) Command.fh Command.tcc \
            Command.xh InterfaceBase.fh InterfaceBase.xh \
            InterfacedBase.xh ParVector.fh \
            ParVector.tcc ParVector.xh Parameter.fh \
            Parameter.tcc Parameter.xh RefVector.fh \
            RefVector.tcc RefVector.xh Reference.fh \
            Reference.tcc Reference.xh Switch.fh \
            Switch.tcc Switch.xh ClassDocumentation.fh Deleted.h

noinst_LTLIBRARIES = libThePEGInter.la

libThePEGInter_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

