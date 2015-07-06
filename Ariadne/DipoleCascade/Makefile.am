mySOURCES = CascadeBase.cc DipoleState.cc Emitter.cc Parton.cc Dipole.cc \
            String.cc CascadeHandler.cc HardSubSys.cc EMDipole.cc \
            RemnantParton.cc ExtendedDipole.cc MECorrBase.cc DYCorr.cc \
	    DISDipole.cc HardRemnant.cc SoftRemnant.cc ReweightBase.cc

DOCFILES = CascadeBase.h DipoleState.h Emitter.h Parton.h Dipole.h String.h \
           CascadeHandler.h HardSubSys.h PartonTraits.h EMDipole.h \
           RemnantParton.h ExtendedDipole.h MECorrBase.h DYCorr.h \
	   DISDipole.h HardRemnant.h SoftRemnant.h ReweightBase.h


INCLUDEFILES = $(DOCFILES) CascadeBase.fh CascadeBase.icc \
               DipoleState.fh DipoleState.icc Emitter.fh Emitter.icc \
               Parton.fh Parton.icc Dipole.fh Dipole.icc String.fh String.icc \
               CascadeHandler.fh CascadeHandler.icc HardSubSys.icc \
               EMDipole.fh EMDipole.icc RemnantParton.fh RemnantParton.icc \
               ExtendedDipole.fh ExtendedDipole.icc MECorrBase.fh \
	       ReweightBase.fh DYCorr.fh DYCorr.icc DISDipole.fh DISDipole.icc \
	       HardRemnant.fh HardRemnant.icc SoftRemnant.fh SoftRemnant.icc

pkglib_LTLIBRARIES = libArCascade.la

# Version info should be updated if any interface or persistent I/O
# function is changed
libArCascade_la_LDFLAGS = -module -version-info 1:0:0
libArCascade_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude
