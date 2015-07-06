// -*- C++ -*-
//
// HepMC3File.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HepMC3File class.
//

#include "HepMC3File.h"
#include <config.h>
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "HepMC/IO/IO_GenEvent.h"

using namespace ThePEG;

HepMC3File::HepMC3File() 
  : _eventNumber(1), _filename(), _unitchoice(),
    _geneventPrecision(16) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
HepMC3File::HepMC3File(const HepMC3File & x) 
  : AnalysisHandler(x), 
    _eventNumber(x._eventNumber),
    _filename(x._filename), _hepmcio(),
    _unitchoice(x._unitchoice), 
    _geneventPrecision(x._geneventPrecision) {}

IBPtr HepMC3File::clone() const {
  return new_ptr(*this);
}

IBPtr HepMC3File::fullclone() const {
  return new_ptr(*this);
}

void HepMC3File::doinitrun() {
  AnalysisHandler::doinitrun();

  // set default filename unless user-specified name exists
  if ( _filename.empty() )
    _filename = generator()->filename() + ".hepmc";

  HepMC::IO_GenEvent * tmpio 
    = new HepMC::IO_GenEvent(_filename.c_str(), ios::out);
  tmpio->set_precision(_geneventPrecision);
  _hepmcio = tmpio;

}

void HepMC3File::dofinish() {
  if (_hepmcio) {
    delete _hepmcio;
    _hepmcio = 0;
  }
  AnalysisHandler::dofinish();
  cout << "\nHepMC3File: generated HepMC output.\n";
}

void HepMC3File::analyze(tEventPtr event, long, int, int) {
  if (event->number() > _eventNumber) return;


  HepMC::GenEvent * hepmc 
    = HepMCConverter<HepMC::GenEvent>::convert(*event, false,
					       GeV, millimeter);
  if (_hepmcio)
    _hepmcio->write_event(*hepmc);
  delete hepmc;
}

void HepMC3File::persistentOutput(PersistentOStream & os) const {
  os << _eventNumber << _filename 
     << _unitchoice << _geneventPrecision;
}

void HepMC3File::persistentInput(PersistentIStream & is, int) {
  is >> _eventNumber >> _filename 
     >> _unitchoice >> _geneventPrecision;
}


ClassDescription<HepMC3File> HepMC3File::initHepMC3File;
// Definition of the static class description member.

void HepMC3File::Init() {

  static ClassDocumentation<HepMC3File> documentation
    ("This analysis handler will output the event record in HepMC format.");

  static Parameter<HepMC3File,long> interfacePrintEvent
    ("PrintEvent",
     "The number of events that should be printed.",
     &HepMC3File::_eventNumber, 1, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<HepMC3File,string> interfaceFilename
    ("Filename", "Name of the output file",
     &HepMC3File::_filename, "");

  static Parameter<HepMC3File,unsigned int> interfacePrecision
    ("Precision",
     "Choice of output precision for the GenEvent format "
     " (as number of digits).",
     &HepMC3File::_geneventPrecision, 16, 6, 16,
     false, false, Interface::limited);
  
  static Switch<HepMC3File,int> interfaceUnits
    ("Units",
     "Unit choice for energy and length",
     &HepMC3File::_unitchoice, 0, false, false);
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
