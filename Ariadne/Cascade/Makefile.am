SUBDIRS = Models

mySOURCES = CascadeBase.cc DipoleBase.cc Parton.cc QCDDipole.cc ColourIndex.cc \
            EMDipole.cc AriadneHandler.cc ReweightBase.cc \
            Emission.cc EmitterBase.cc DipoleState.cc \
            EmissionGenerator.cc RemnantParton.cc \
            History.cc BornCheckerBase.cc ScaleSetter.cc \
            ResonanceFinder.cc Junction.cc ResonanceParton.cc \
            QCDDipoleFinder.cc EMDipoleFinder.cc Resonance.cc \
            DISFinder.cc ConsistencyChecker.cc StateDipole.cc

DOCFILES = CascadeBase.h DipoleBase.h Parton.h QCDDipole.h ColourIndex.h \
           EMDipole.h AriadneHandler.h ReweightBase.h Emission.h \
           EmitterBase.h DipoleState.h PartonTraits.h \
           EmissionGenerator.h RemnantParton.h \
           History.h BornCheckerBase.h ScaleSetter.h ResonanceFinder.h \
           Junction.h ResonanceParton.h QCDDipoleFinder.h \
           EMDipoleFinder.h Resonance.h DISFinder.h ConsistencyChecker.h \
           StateDipole.h

INCLUDEFILES = $(DOCFILES) CascadeBase.fh DipoleState.fh AriadneHandler.fh \
               DipoleBase.fh Parton.fh QCDDipole.fh Emission.fh \
               EMDipole.fh ReweightBase.fh EmitterBase.fh RemnantParton.fh \
               History.fh BornCheckerBase.fh \
               ScaleSetter.fh ResonanceFinder.fh Junction.fh \
               ResonanceParton.fh QCDDipoleFinder.fh \
               EMDipoleFinder.fh Resonance.fh DISFinder.fh ConsistencyChecker.fh \
               StateDipole.fh

pkglib_LTLIBRARIES = libAriadne5.la MinPTScale.la MaxPTScale.la

# Version info should be updated if any interface or persistent I/O
# function is changed
libAriadne5_la_LDFLAGS = -module -version-info 1:0:0
libAriadne5_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
libAriadne5_la_LIBADD = Models/libAriadne5Models.la ../Config/libAriadne5Config.la

MinPTScale_la_SOURCES = MinPTScale.cc MinPTScale.h
MinPTScale_la_LDFLAGS = -module -version-info 1:0:0

MaxPTScale_la_SOURCES = MaxPTScale.cc MaxPTScale.h
MaxPTScale_la_LDFLAGS = -module -version-info 1:0:0

include $(top_srcdir)/Config/Makefile.aminclude

..check:
	cd ..; make check
