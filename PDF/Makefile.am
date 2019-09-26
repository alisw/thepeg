mySOURCES = PartonExtractor.cc PDFBase.cc RemnantHandler.cc \
            BeamParticleData.cc PolarizedBeamParticleData.cc \
            PartonBin.cc NoPDF.cc NoRemnants.cc PDFCuts.cc \
            PartonBinInstance.cc SoftRemnantHandler.cc

DOCFILES = BeamParticleData.h PolarizedBeamParticleData.h \
           NoPDF.h NoRemnants.h PDFBase.h PDFCuts.h \
           PartonBin.h PartonExtractor.h RemnantHandler.h \
           PDF.h PartonBinInstance.h SoftRemnantHandler.h

INCLUDEFILES = $(DOCFILES) BeamParticleData.xh \
               PDFBase.xh PartonExtractor.fh \
               PartonExtractor.xh RemnantHandler.xh \
               PDF.fh PartonBin.fh PolarizedBeamParticleData.fh

INPUTFILES = PDFsets.index

CLEANFILES = .done-all-links

noinst_LTLIBRARIES = libThePEGPDF.la
pkglib_LTLIBRARIES = LeptonLeptonPDF.la UnResolvedRemnant.la \
          GRVBase.la GRV94L.la GRV94M.la WeizsackerWilliamsPDF.la \
          BudnevPDF.la

if USELHAPDF
pkglib_LTLIBRARIES += ThePEGLHAPDF.la
LHAPDF_SOURCES = LHAPDF6.cc LHAPDF6.h

ThePEGLHAPDF_la_LDFLAGS =  $(AM_LDFLAGS) $(LHAPDF_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
ThePEGLHAPDF_la_LIBADD = $(LHAPDF_LIBS)
ThePEGLHAPDF_la_SOURCES = $(LHAPDF_SOURCES)
ThePEGLHAPDF_la_CPPFLAGS = $(AM_CPPFLAGS) $(LHAPDF_CPPFLAGS)
endif

libThePEGPDF_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
libThePEGPDF_la_CPPFLAGS = $(AM_CPPFLAGS)

# Version info should be updated if any interface or persistent I/O
# function is changed
LeptonLeptonPDF_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
LeptonLeptonPDF_la_SOURCES = LeptonLeptonPDF.cc LeptonLeptonPDF.h

# Version info should be updated if any interface or persistent I/O
# function is changed
WeizsackerWilliamsPDF_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
WeizsackerWilliamsPDF_la_SOURCES = WeizsackerWilliamsPDF.cc WeizsackerWilliamsPDF.h

# Version info should be updated if any interface or persistent I/O
# function is changed
BudnevPDF_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
BudnevPDF_la_SOURCES = BudnevPDF.cc BudnevPDF.h


# Version info should be updated if any interface or persistent I/O
# function is changed
UnResolvedRemnant_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
UnResolvedRemnant_la_SOURCES = UnResolvedRemnant.cc UnResolvedRemnant.h

# Version info should be updated if any interface or persistent I/O
# function is changed
GRVBase_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
GRVBase_la_SOURCES = GRVBase.cc GRVBase.h

# Version info should be updated if any interface or persistent I/O
# function is changed
GRV94L_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
GRV94L_la_SOURCES = GRV94L.cc GRV94L.h

# Version info should be updated if any interface or persistent I/O
# function is changed
GRV94M_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
GRV94M_la_SOURCES = GRV94M.cc GRV94M.h

dist_pkgdata_DATA = PDFsets.index

all-local: .done-all-links

.done-all-links:
@EMPTY@ifdef SHOWCOMMAND
	for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	echo "stamp" > .done-all-links
@EMPTY@else
	@echo "sym-linking input files files..."
	@for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	@echo "stamp" > .done-all-links
@EMPTY@endif

include $(top_srcdir)/Config/Makefile.aminclude
