mySOURCES = Event.cc Collision.cc SubProcess.cc SubProcessGroup.cc Step.cc Particle.cc \
          EventInfoBase.cc ColourLine.cc ColourBase.cc SpinInfo.cc \
          EventConfig.cc ColourSinglet.cc RemnantParticle.cc MultiColour.cc \
	  HelicityVertex.cc

DOCFILES = EventConfig.h Collision.h Event.h Particle.h ParticleTraits.h \
           SelectorBase.h StandardSelectors.h Step.h SubProcess.h SubProcessGroup.h \
           EventInfoBase.h ColourLine.h ColourBase.h SpinInfo.h \
           ColourSinglet.h TmpTransform.h RemnantParticle.h MultiColour.h \
	   HelicityVertex.h RhoDMatrix.h

INCLUDEFILES = $(DOCFILES) Collision.tcc \
               Particle.fh Particle.tcc \
               Step.tcc SubProcess.tcc RemnantParticle.fh HelicityVertex.fh

noinst_LTLIBRARIES = libThePEGEvent.la
libThePEGEvent_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

