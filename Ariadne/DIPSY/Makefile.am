AUTOMAKE_OPTIONS = -Wno-portability

mySOURCES = DipoleEventHandler.cc WaveFunction.cc SimpleProton.cc \
            VirtualPhoton.cc DipoleState.cc Dipole.cc Parton.cc DipoleXSec.cc \
            ImpactParameterGenerator.cc WFInfo.cc PhotonWFInfo.cc \
            VectorMesonBase.cc PhotonDipoleState.cc SimpleProtonState.cc \
            Emitter.cc Swinger.cc  DipoleAnalysisHandler.cc \
            TotalXSecAnalysis.cc EventFiller.cc DipoleAbsorber.cc \
            SmallDipoleAbsorber.cc EffectiveParton.cc RealParton.cc \
            RealPartonState.cc DiffractiveEventFiller.cc ParticleInfo.cc \
            NucleusData.cc

DOCFILES = DipoleEventHandler.h WaveFunction.h SimpleProton.h VirtualPhoton.h \
           DipoleState.h Dipole.h Parton.h DipoleXSec.h ImpactParameters.h \
           ImpactParameterGenerator.h WFInfo.h PhotonWFInfo.h VectorMesonBase.h \
           PhotonDipoleState.h SimpleProtonState.h Emitter.h Swinger.h \
           DipoleAnalysisHandler.h TotalXSecAnalysis.h EventFiller.h \
           DipoleAbsorber.h SmallDipoleAbsorber.h EffectiveParton.h \
           RealParton.h RealPartonState.h DiffractiveEventFiller.h ParticleInfo.h \
           NucleusData.h

INCLUDEFILES = $(DOCFILES) DipoleEventHandler.fh \
               WaveFunction.icc WaveFunction.fh \
               VirtualPhoton.icc DipoleState.fh DipoleState.icc Dipole.fh \
               Dipole.icc Parton.fh Parton.icc DipoleXSec.fh DipoleXSec.icc \
               ImpactParameters.icc ImpactParameterGenerator.fh  WFInfo.fh \
               PhotonWFInfo.icc VectorMesonBase.icc PhotonDipoleState.icc \
               SimpleProtonState.icc Emitter.fh Emitter.icc Swinger.fh \
               Swinger.icc DipoleAnalysisHandler.fh EventFiller.fh \
               DipoleAbsorber.fh EffectiveParton.fh \
               RealParton.fh RealPartonState.fh CPUTimer.h

pkglib_LTLIBRARIES = libDIPSY.la OldStyleEmitter.la PT1DEmitter.la \
                     ElasticXSecAnalysis.la LargePTDipoleAbsorber.la \
                     SimpleNucleus.la RecoilSwinger.la \
                     FSDipole5Ordering.la FSAnalysis.la GapAnalysis.la HIAnalysis.la \
                     FixedImpactGenerator.la NoFluctXSecAnalysis.la GlauberAnalysis.la \
                     AnalysisProgress.la SemiInclusiveXSecAnalysis.la \
                     GaussianImpactGenerator.la DipoleDensityAnalysis.la

INPUTFILES = TestXSecs.in RivetAnalyses.in DIPSYRemove.in DIPSYDefaults.in CurrentTune.in Tune27.in

CLEANFILES = done-all-links

dist_pkgdata_DATA = $(INPUTFILES)
rpodir = $(pkglibdir)
nodist_rpo_DATA = DIPSYDefaults.rpo

# Version info should be updated if any interface or persistent I/O
# function is changed
libDIPSY_la_LDFLAGS = -module -version-info 1:0:0
# libDIPSY_la_LIBADD = ../Config/libAriadne5Config.la
libDIPSY_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

OldStyleEmitter_la_LDFLAGS = -module -version-info 1:0:0
OldStyleEmitter_la_SOURCES = OldStyleEmitter.cc OldStyleEmitter.h

FSDipole5Ordering_la_LDFLAGS = -module -version-info 1:0:0
FSDipole5Ordering_la_SOURCES = FSDipole5Ordering.cc FSDipole5Ordering.h

FSAnalysis_la_LDFLAGS = -module -version-info 1:0:0
FSAnalysis_la_SOURCES = FSAnalysis.cc FSAnalysis.h

GapAnalysis_la_LDFLAGS = -module -version-info 1:0:0
GapAnalysis_la_SOURCES = GapAnalysis.cc GapAnalysis.h

HIAnalysis_la_LDFLAGS = -module -version-info 1:0:0
HIAnalysis_la_SOURCES = HIAnalysis.cc HIAnalysis.h

PT1DEmitter_la_LDFLAGS = -module -version-info 1:0:0
PT1DEmitter_la_SOURCES = PT1DEmitter.cc PT1DEmitter.h

ElasticXSecAnalysis_la_LDFLAGS = -module -version-info 1:0:0
ElasticXSecAnalysis_la_SOURCES = ElasticXSecAnalysis.cc ElasticXSecAnalysis.h

DipoleDensityAnalysis_la_LDFLAGS = -module -version-info 1:0:0
DipoleDensityAnalysis_la_SOURCES = DipoleDensityAnalysis.cc DipoleDensityAnalysis.h

NoFluctXSecAnalysis_la_LDFLAGS = -module -version-info 1:0:0
NoFluctXSecAnalysis_la_SOURCES = NoFluctXSecAnalysis.cc NoFluctXSecAnalysis.h

SemiInclusiveXSecAnalysis_la_LDFLAGS = -module -version-info 1:0:0
SemiInclusiveXSecAnalysis_la_SOURCES = SemiInclusiveXSecAnalysis.cc SemiInclusiveXSecAnalysis.h

GlauberAnalysis_la_LDFLAGS = -module -version-info 1:0:0
GlauberAnalysis_la_SOURCES = GlauberAnalysis.cc GlauberAnalysis.h

AnalysisProgress_la_LDFLAGS = -module -version-info 1:0:0
AnalysisProgress_la_SOURCES = AnalysisProgress.cc AnalysisProgress.h

LargePTDipoleAbsorber_la_LDFLAGS = -module -version-info 1:0:0
LargePTDipoleAbsorber_la_SOURCES = LargePTDipoleAbsorber.cc LargePTDipoleAbsorber.h

SimpleNucleus_la_LDFLAGS = -module -version-info 1:0:0
SimpleNucleus_la_SOURCES = SimpleNucleus.cc SimpleNucleus.h \
                           SimpleNucleusState.cc SimpleNucleusState.h

RecoilSwinger_la_LDFLAGS = -module -version-info 1:0:0
RecoilSwinger_la_SOURCES = RecoilSwinger.cc RecoilSwinger.h

FixedImpactGenerator_la_LDFLAGS = -module -version-info 1:0:0
FixedImpactGenerator_la_SOURCES = FixedImpactGenerator.cc FixedImpactGenerator.h

GaussianImpactGenerator_la_LDFLAGS = -module -version-info 1:0:0
GaussianImpactGenerator_la_SOURCES = GaussianImpactGenerator.cc GaussianImpactGenerator.h

EXTRACLASSES = OldStyleEmitter.la PT1DEmitter.la ElasticXSecAnalysis.la \
               LargePTDipoleAbsorber.la RecoilSwinger.la SimpleNucleus.la \
               FSAnalysis.la GapAnalysis.la HIAnalysis.la \
               FixedImpactGenerator.la FSDipole5Ordering.la

done-all-links:
@EMPTY@ifdef SHOWCOMMAND
	for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	echo "stamp" > done-all-links
@EMPTY@else
	@echo "sym-linking input files files..."
	@for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	@echo "stamp" > done-all-links
@EMPTY@endif

DIPSYDefaults.rpo: done-all-links DIPSYDefaults.in libDIPSY.la $(EXTRACLASSES) $(THEPEGLIB)/ThePEGDefaults.rpo ../lib/Ariadne5Defaults.rpo
	$(SETUPTHEPEG) -L../lib -L$(THEPEGLIB) -L ../../TheP8I/lib -L .libs -L ../src/.libs --exitonerror --init -r ../lib/Ariadne5Defaults.rpo -o DIPSYDefaults.rpo DIPSYDefaults.in
	cd ../lib; $(LN_S) -f ../DIPSY/DIPSYDefaults.rpo .


Test5%.run: Test5%.in DIPSYDefaults.rpo CurrentTune.in RivetAnalyses.in
	$(SETUPTHEPEG) --exitonerror -L .libs -L ../../Pythia7/src/.libs -r DIPSYDefaults.rpo $<

%.run: %.in DIPSYDefaults.rpo CurrentTune.in RivetAnalyses.in
	$(SETUPTHEPEG) -L../lib --exitonerror -L .libs -r DIPSYDefaults.rpo $<

%.out: %.run
	time $(RUNTHEPEG) --tics -d 0 $<

valgrind:
	valgrind --leak-check=full --num-callers=25 --track-fds=yes --freelist-vol=100000000 --leak-resolution=med --trace-children=yes $(top_builddir)/../ThePEG/src/setupThePEG -L ../../ThePEG/lib -L../../Pythia7/lib --exitonerror -L../../TheP8I/lib -L .libs -L ../lib -r ../lib/Ariadne5Defaults.rpo TestFull.in &> /tmp/valgrind.out
	valgrind --leak-check=full --num-callers=25 --track-fds=yes --freelist-vol=100000000 --leak-resolution=med --trace-children=yes $(top_builddir)/../ThePEG/src/runThePEG -N 10 TestFull.run >> /tmp/valgrind.out 2>&1

install-data-local:
	LD_LIBRARY_PATH=$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH $(DESTDIR)$(bindir)/setupThePEG --exitonerror --init -i $(DESTDIR)$(pkgdatadir) -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(srcdir)/DIPSYDefaults.in

unregister-from-repo:
	LD_LIBRARY_PATH=$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH $(DESTDIR)$(bindir)/setupThePEG --exitonerror --init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(srcdir)/DIPSYRemove.in



include $(top_srcdir)/Config/Makefile.aminclude
