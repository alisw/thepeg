mySOURCES = LorentzRotation.cc SpinHalfLorentzRotation.cc SpinOneLorentzRotation.cc

DOCFILES = Lorentz5Vector.h LorentzRotation.h LorentzVector.h \
           ThreeVector.h Transverse.h SpinOneLorentzRotation.h \
           SpinHalfLorentzRotation.h HepMCTraits.h HepMCConverter.h

INCLUDEFILES = $(DOCFILES) Lorentz5Vector.fh \
               LorentzRotation.fh LorentzVector.fh \
               ThreeVector.fh Transverse.fh \
               SpinOneLorentzRotation.fh \
               HepMCConverter.tcc

noinst_LTLIBRARIES = libThePEGVectors.la

libThePEGVectors_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

