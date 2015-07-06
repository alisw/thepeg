mySOURCES = FSGluonEmission.cc FSGluonEmitter.cc FSQQEmission.cc FSQQEmitter.cc \
            ColourChargeRegions.cc SpecialGluonEmission.cc \
            SpecialGluonEmitter.cc RemnantModel.cc ColourResonanceModel.cc \
            RemnantGluonEmission.cc ISGtoQEmission.cc ISQEmission.cc \
            ISQEmitter.cc ISGtoQEmitter.cc ISQtoGEmitter.cc ISQtoGEmission.cc \
            DISME.cc RemnantGluonEmitter.cc SpecialQQEmission.cc \
            SpecialQQEmitter.cc RRGluonEmitter.cc RRGluonEmission.cc DYME.cc \
            RecoilGluonEmitter.cc DipoleSwinger.cc DipoleSwing.cc

DOCFILES = FSGluonEmission.h FSGluonEmitter.h FSQQEmission.h FSQQEmitter.h \
           ColourChargeRegions.h SpecialGluonEmission.h SpecialGluonEmitter.h \
           RemnantModel.h PseudoParton.h ColourResonanceModel.h \
           RemnantGluonEmission.h ISGtoQEmission.h ISQEmission.h ISQEmitter.h \
           ISGtoQEmitter.h ISQtoGEmitter.h ISQtoGEmission.h DISME.h \
           RemnantGluonEmitter.h SpecialQQEmission.h SpecialQQEmitter.h DYME.h \
           RRGluonEmitter.h RRGluonEmission.h RecoilGluonEmitter.h \
           DipoleSwinger.h DipoleSwing.h

INCLUDEFILES = $(DOCFILES) RemnantModel.fh ColourResonanceModel.fh

noinst_LTLIBRARIES = libAriadne5Models.la

libAriadne5Models_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)




include $(top_srcdir)/Config/Makefile.aminclude

..debug:
	cd ..; make debug

..check:
	cd ../..; make check

