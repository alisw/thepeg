AIDAWRAPPERS = LWH/AIAnalysisFactory.h LWH/AIAxis.h LWH/AIHistogram1D.h LWH/AIHistogram2D.h \
               LWH/AIHistogramFactory.h LWH/AITree.h LWH/AITreeFactory.h \
               LWH/AIManagedObject.h LWH/AIDataPoint.h LWH/AIDataPointSet.h \
               LWH/AIDataPointSetFactory.h LWH/AIMeasurement.h

LWHHEADERS = LWH/AnalysisFactory.h LWH/Axis.h LWH/Tree.h LWH/TreeFactory.h \
             LWH/Histogram1D.h LWH/Histogram2D.h LWH/HistogramFactory.h LWH/ManagedObject.h \
             LWH/VariAxis.h LWH/DataPoint.h LWH/DataPointSet.h LWH/DataPointSetFactory.h \
             LWH/Measurement.h

mySOURCES = LWHFactory.cc

pkglib_LTLIBRARIES = LWHFactory.la XSecCheck.la ProgressLog.la

if HAVE_RIVET
  DOCFILES = RivetAnalysis.h NLORivetAnalysis.h
  pkglib_LTLIBRARIES += RivetAnalysis.la
  RivetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(RIVETINCLUDE) $(HEPMCINCLUDE)
  RivetAnalysis_la_LIBADD = $(HEPMCLIBS) $(RIVETLIBS)
  RivetAnalysis_la_SOURCES = RivetAnalysis.cc NLORivetAnalysis.cc $(INCLUDEFILES)
  RivetAnalysis_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
else
  DOCFILES = $(AIDAWRAPPERS) $(LWHHEADERS)
endif
DOCFILES += LWHFactory.h

INCLUDEFILES = $(DOCFILES) 

LWHINCLUDE = -I$(top_builddir)/include/ThePEG/Analysis/LWH

LWHFactory_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
LWHFactory_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
LWHFactory_la_CPPFLAGS = $(AM_CPPFLAGS) $(LWHINCLUDE)

XSecCheck_la_SOURCES = XSecCheck.cc XSecCheck.h
XSecCheck_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)

ProgressLog_la_SOURCES = ProgressLog.cc ProgressLog.h
ProgressLog_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)

noinst_LTLIBRARIES = libThePEGHist.la
libThePEGHist_la_SOURCES = \
  FactoryBase.cc FactoryBase.fh FactoryBase.h AIDA_helper.h
libThePEGHist_la_CPPFLAGS = $(AM_CPPFLAGS) $(LWHINCLUDE)

if HAVE_HEPMC
  pkglib_LTLIBRARIES += HepMCAnalysis.la
  HepMCAnalysis_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
  HepMCAnalysis_la_SOURCES = HepMCFile.h    HepMCFile.cc    NLOHepMCFile.h    NLOHepMCFile.cc	HIHepMCFile.h HIHepMCFile.cc
  HepMCAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(HEPMCINCLUDE)
  HepMCAnalysis_la_LIBADD = $(HEPMCLIBS)
endif

pkglib_LTLIBRARIES += GraphvizPlot.la
GraphvizPlot_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
GraphvizPlot_la_SOURCES = GraphvizPlot.h GraphvizPlot.cc

LWH.tgz: $(AIDAWRAPPERS) $(LWHHEADERS)
	tar czf LWH.tgz  $(AIDAWRAPPERS) $(LWHHEADERS)

include $(top_srcdir)/Config/Makefile.aminclude






