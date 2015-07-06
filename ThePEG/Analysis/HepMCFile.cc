// -*- C++ -*-
//
// HepMCFile.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HepMCFile class.
//

#include "HepMCFile.h"
#include <config.h>
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"

using namespace ThePEG;

HepMCFile::HepMCFile() 
  : _eventNumber(1), _format(1), _filename(), _unitchoice(),
    _geneventPrecision(16) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
HepMCFile::HepMCFile(const HepMCFile & x) 
  : AnalysisHandler(x), 
    _eventNumber(x._eventNumber), _format(x._format), 
    _filename(x._filename), _hepmcio(), _hepmcdump(), 
    _unitchoice(x._unitchoice), 
    _geneventPrecision(x._geneventPrecision) {}

IBPtr HepMCFile::clone() const {
  return new_ptr(*this);
}

IBPtr HepMCFile::fullclone() const {
  return new_ptr(*this);
}

void HepMCFile::doinitrun() {
  AnalysisHandler::doinitrun();

  // set default filename unless user-specified name exists
  if ( _filename.empty() )
    _filename = generator()->filename() + ".hepmc";

  switch ( _format ) {
  default: {
    HepMC::IO_GenEvent * tmpio 
      = new HepMC::IO_GenEvent(_filename.c_str(), ios::out);
    tmpio->precision(_geneventPrecision);
    _hepmcio = tmpio;
    break;
  }
  case 2: 
    _hepmcio = new HepMC::IO_AsciiParticles(_filename.c_str(), ios::out); 
    break;
  case 5: 
    _hepmcio = 0; 
    _hepmcdump.open(_filename.c_str()); 
    break;
  }
}

void HepMCFile::dofinish() {
  if (_hepmcio) {
    delete _hepmcio;
    _hepmcio = 0;
  }
  else
    _hepmcdump.close();
  AnalysisHandler::dofinish();
  cout << "\nHepMCFile: generated HepMC output.\n";
}

void HepMCFile::analyze(tEventPtr event, long, int, int) {
  if (event->number() > _eventNumber) return;

  Energy eUnit;
  Length lUnit;
  switch (_unitchoice) {
  default: eUnit = GeV; lUnit = millimeter; break;
  case 1:  eUnit = MeV; lUnit = millimeter; break;
  case 2:  eUnit = GeV; lUnit = centimeter; break;
  case 3:  eUnit = MeV; lUnit = centimeter; break;
  }

  HepMC::GenEvent * hepmc 
    = HepMCConverter<HepMC::GenEvent>::convert(*event, false,
					       eUnit, lUnit);
  if (_hepmcio)
    _hepmcio->write_event(hepmc);
  else
    hepmc->print(_hepmcdump);
  delete hepmc;
}

void HepMCFile::persistentOutput(PersistentOStream & os) const {
  os << _eventNumber << _format << _filename 
     << _unitchoice << _geneventPrecision;
}

void HepMCFile::persistentInput(PersistentIStream & is, int) {
  is >> _eventNumber >> _format >> _filename 
     >> _unitchoice >> _geneventPrecision;
}


ClassDescription<HepMCFile> HepMCFile::initHepMCFile;
// Definition of the static class description member.

void HepMCFile::Init() {

  static ClassDocumentation<HepMCFile> documentation
    ("This analysis handler will output the event record in HepMC format.");

  static Parameter<HepMCFile,long> interfacePrintEvent
    ("PrintEvent",
     "The number of events that should be printed.",
     &HepMCFile::_eventNumber, 1, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<HepMCFile,int> interfaceFormat
    ("Format",
     "Output format (1 = GenEvent, 2 = AsciiParticles, 5 = HepMC dump)",
     &HepMCFile::_format, 1, false, false);
  static SwitchOption interfaceFormatGenEvent
    (interfaceFormat,
     "GenEvent",
     "IO_GenEvent format",
     1);
  static SwitchOption interfaceFormatAsciiParticles
    (interfaceFormat,
     "AsciiParticles",
     "Deprecated (IO_AsciiParticles format)",
     2);
  static SwitchOption interfaceFormatDump
    (interfaceFormat,
     "Dump",
     "Event dump (human readable)",
     5);

  static Parameter<HepMCFile,string> interfaceFilename
    ("Filename", "Name of the output file",
     &HepMCFile::_filename, "");

  static Parameter<HepMCFile,unsigned int> interfacePrecision
    ("Precision",
     "Choice of output precision for the GenEvent format "
     " (as number of digits).",
     &HepMCFile::_geneventPrecision, 16, 6, 16,
     false, false, Interface::limited);
  
  static Switch<HepMCFile,int> interfaceUnits
    ("Units",
     "Unit choice for energy and length",
     &HepMCFile::_unitchoice, 0, false, false);
  static SwitchOption interfaceUnitsGeV_mm
    (interfaceUnits,
     "GeV_mm",
     "Use GeV and mm as units.",
     0);
  static SwitchOption interfaceUnitsMeV_mm
    (interfaceUnits,
     "MeV_mm",
     "Use MeV and mm as units.",
     1);
  static SwitchOption interfaceUnitsGeV_cm
    (interfaceUnits,
     "GeV_cm",
     "Use GeV and cm as units.",
     2);
  static SwitchOption interfaceUnitsMeV_cm
    (interfaceUnits,
     "MeV_cm",
     "Use MeV and cm as units.",
     3);
}
