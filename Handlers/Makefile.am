mySOURCES = EventHandler.cc SubProcessHandler.cc \
            HandlerGroup.cc Hint.cc XComb.cc AnalysisHandler.cc \
            CascadeHandler.cc MultipleInteractionHandler.cc DecayHandler.cc \
            EventManipulator.cc StepHandler.cc HadronizationHandler.cc \
            HandlerBase.cc \
            FlavourGenerator.cc LuminosityFunction.cc \
            StandardEventHandler.cc  \
            SamplerBase.cc  ClusterCollapser.cc \
            PtGenerator.cc ZGenerator.cc   \
            StandardXComb.cc StdXCombGroup.cc

DOCFILES = AnalysisHandler.h CascadeHandler.h \
           DecayHandler.h EventHandler.h EventManipulator.h \
           FlavourGenerator.h HadronizationHandler.h \
           HandlerBase.h HandlerGroup.h Hint.h \
           LuminosityFunction.h \
           MultipleInteractionHandler.h StandardEventHandler.h \
           StepHandler.h SubProcessHandler.h XComb.h LastXCombInfo.h \
           SamplerBase.h  ClusterCollapser.h \
           PtGenerator.h ZGenerator.h StandardXComb.h \
	   StdXCombGroup.h

INCLUDEFILES = $(DOCFILES) HandlerGroup.tcc \
               StandardEventHandler.fh \
               SubProcessHandler.fh EventHandler.fh SamplerBase.fh \
               ClusterCollapser.fh StandardXComb.fh \
	       StdXCombGroup.fh AnalysisHandler.fh


noinst_LTLIBRARIES = libThePEGHandlers.la
pkglib_LTLIBRARIES = FixedCMSLuminosity.la FixedTargetLuminosity.la \
          ACDCSampler.la SimpleFlavour.la GaussianPtGenerator.la \
          SimpleZGenerator.la


libThePEGHandlers_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

# Version info should be updated if any interface or persistent I/O
# function is changed
FixedCMSLuminosity_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
FixedCMSLuminosity_la_SOURCES = FixedCMSLuminosity.cc FixedCMSLuminosity.h

# Version info should be updated if any interface or persistent I/O
# function is changed
FixedTargetLuminosity_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
FixedTargetLuminosity_la_SOURCES = FixedTargetLuminosity.cc FixedTargetLuminosity.h

# Version info should be updated if any interface or persistent I/O
# function is changed
ACDCSampler_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
ACDCSampler_la_SOURCES = ACDCSampler.cc ACDCSampler.h

# Version info should be updated if any interface or persistent I/O
# function is changed
SimpleFlavour_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
SimpleFlavour_la_SOURCES = SimpleFlavour.cc SimpleFlavour.h

# Version info should be updated if any interface or persistent I/O
# function is changed
GaussianPtGenerator_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
GaussianPtGenerator_la_SOURCES = GaussianPtGenerator.cc GaussianPtGenerator.h

# Version info should be updated if any interface or persistent I/O
# function is changed
SimpleZGenerator_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
SimpleZGenerator_la_SOURCES = SimpleZGenerator.cc SimpleZGenerator.h

include $(top_srcdir)/Config/Makefile.aminclude

