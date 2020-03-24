SUBDIRS = WaveFunction Vertex

mySOURCES = LorentzSpinor.cc \
          LorentzSpinorBar.cc LorentzTensor.cc \
          FermionSpinInfo.cc VectorSpinInfo.cc ScalarSpinInfo.cc TensorSpinInfo.cc \
          LorentzRSSpinor.cc LorentzRSSpinorBar.cc RSFermionSpinInfo.cc

DOCFILES = LorentzPolarizationVector.h LorentzSpinor.h \
           LorentzSpinorBar.h LorentzTensor.h \
           FermionSpinInfo.h VectorSpinInfo.h ScalarSpinInfo.h TensorSpinInfo.h \
           LorentzRSSpinor.h LorentzRSSpinorBar.h RSFermionSpinInfo.h \
	   epsilon.h HelicityFunctions.h

INCLUDEFILES = $(DOCFILES) \
               LorentzSpinor.fh LorentzSpinor.tcc LorentzSpinorBar.tcc \
               LorentzSpinorBar.fh \
               FermionSpinInfo.fh \
               LorentzTensor.fh \
               VectorSpinInfo.fh \
               ScalarSpinInfo.fh \
               TensorSpinInfo.fh LorentzTensor.tcc \
               HelicityDefinitions.h \
               LorentzRSSpinor.fh LorentzRSSpinorBar.fh RSFermionSpinInfo.fh \
               LorentzRSSpinor.tcc LorentzRSSpinorBar.tcc

noinst_LTLIBRARIES = libThePEGHelicity.la

libThePEGHelicity_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
libThePEGHelicity_la_LIBADD = Vertex/libThePEGVertex.la \
			      WaveFunction/libThePEGWaveFunction.la

include $(top_srcdir)/Config/Makefile.aminclude

