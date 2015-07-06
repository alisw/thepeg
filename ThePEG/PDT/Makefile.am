mySOURCES = ParticleData.cc DecayMode.cc MatcherBase.cc Matcher.cc Decayer.cc \
          MassGenerator.cc WidthGenerator.cc DummyDecayer.cc \
          ConstituentParticleData.cc FlatDecayer.cc PDT.cc RemnantData.cc \
          RemnantDecayer.cc SimpleBaryonRemnantDecayer.cc MixedParticleData.cc

DOCFILES = CombinedMatcher.h DecayMode.h WidthGenerator.h Decayer.h \
           EnumParticles.h MassGenerator.h Matcher.h MatcherBase.h PDT.h \
           PID.h ParticleData.h StandardMatchers.h DummyDecayer.h \
           ConstituentParticleData.h FlatDecayer.h RemnantData.h \
           RemnantDecayer.h SimpleBaryonRemnantDecayer.h MixedParticleData.h

INCLUDEFILES = $(DOCFILES) DecayMode.fh \
               DecayMode.xh WidthGenerator.fh Decayer.fh \
               MassGenerator.fh \
               Matcher.tcc ParticleData.xh \
               RemnantData.fh RemnantDecayer.fh MixedParticleData.fh

noinst_LTLIBRARIES = libThePEGPDT.la
pkglib_LTLIBRARIES = BreitWignerMass.la OmegaPhi3PiDecayer.la \
          V2PPDecayer.la DalitzDecayer.la Tau2HadronsDecayer.la \
          Onium3GDecayer.la QuarksToHadronsDecayer.la ColourPairDecayer.la \
          WeakToHadronsDecayer.la

libThePEGPDT_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

# Version info should be updated if any interface or persistent I/O
# function is changed
BreitWignerMass_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
BreitWignerMass_la_SOURCES = BreitWignerMass.cc BreitWignerMass.h

# Version info should be updated if any interface or persistent I/O
# function is changed
OmegaPhi3PiDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
OmegaPhi3PiDecayer_la_SOURCES = OmegaPhi3PiDecayer.cc OmegaPhi3PiDecayer.h

# Version info should be updated if any interface or persistent I/O
# function is changed
V2PPDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
V2PPDecayer_la_SOURCES = V2PPDecayer.cc V2PPDecayer.h

# Version info should be updated if any interface or persistent I/O
# function is changed
DalitzDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
DalitzDecayer_la_SOURCES = DalitzDecayer.cc DalitzDecayer.h

# Version info should be updated if any interface or persistent I/O
# function is changed
Tau2HadronsDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
Tau2HadronsDecayer_la_SOURCES = Tau2HadronsDecayer.cc Tau2HadronsDecayer.h

# Version info should be updated if any interface or persistent I/O
# function is changed
Onium3GDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
Onium3GDecayer_la_SOURCES = Onium3GDecayer.cc Onium3GDecayer.h

# Version info should be updated if any interface or persistent I/O
# function is changed
QuarksToHadronsDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
QuarksToHadronsDecayer_la_SOURCES = QuarksToHadronsDecayer.cc QuarksToHadronsDecayer.h

# Version info should be updated if any interface or persistent I/O
# function is changed
ColourPairDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
ColourPairDecayer_la_SOURCES = ColourPairDecayer.cc ColourPairDecayer.h

# Version info should be updated if any interface or persistent I/O
# function is changed
WeakToHadronsDecayer_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
WeakToHadronsDecayer_la_SOURCES = WeakToHadronsDecayer.cc WeakToHadronsDecayer.h

include $(top_srcdir)/Config/Makefile.aminclude

