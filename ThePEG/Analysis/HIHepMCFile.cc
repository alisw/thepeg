// -*- C++ -*-
//
// HIHepMCFile.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HIHepMCFile class.
//

#include "HIHepMCFile.h"
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
#include "HepMC/GenCrossSection.h"

using namespace ThePEG;

HIHepMCFile::HIHepMCFile() 
  : _eventNumber(1), _format(1), _filename(), _unitchoice(),
    _geneventPrecision(16) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
HIHepMCFile::HIHepMCFile(const HIHepMCFile & x) 
  : AnalysisHandler(x), 
    _eventNumber(x._eventNumber), _format(x._format), 
    _filename(x._filename), _hepmcio(), _hepmcdump(), 
    _unitchoice(x._unitchoice), 
    _geneventPrecision(x._geneventPrecision) {}

IBPtr HIHepMCFile::clone() const {
  return new_ptr(*this);
}

IBPtr HIHepMCFile::fullclone() const {
  return new_ptr(*this);
}

void HIHepMCFile::doinitrun() {
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

void HIHepMCFile::dofinish() {
  if (_hepmcio) {
    delete _hepmcio;
    _hepmcio = 0;
  }
  else
    _hepmcdump.close();
  AnalysisHandler::dofinish();
  cout << "\nHIHepMCFile: generated HepMC output.\n";
}

void HIHepMCFile::analyze(tEventPtr event, long, int, int) {
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

  HepMC::HeavyIon heavyion(1,1,1,1,1,1);
  
  const LorentzPoint v1 = event->incoming().first->vertex();
  const LorentzPoint v2 = event->incoming().second->vertex();
  
  double bpar = (v1 - v2).perp()/femtometer;
  heavyion.HepMC::HeavyIon::set_event_plane_angle(atan2((v1 - v2).y(),(v1 - v2).x()));
  heavyion.HepMC::HeavyIon::set_impact_parameter(float(bpar));  

  // Clear and blatant abuse of the Pdf info container!!
  HepMC::PdfInfo pdfinfo(1,1,event->optionalWeight("averageKappa"),event->optionalWeight("junctions"),event->optionalWeight("lambdaSum"),1,1);

  
  hepmc->set_heavy_ion(heavyion);
  hepmc->set_pdf_info(pdfinfo);
  
  if (_hepmcio)
    _hepmcio->write_event(hepmc);
  else
    hepmc->print(_hepmcdump);
  delete hepmc;
}

void HIHepMCFile::persistentOutput(PersistentOStream & os) const {
  os << _eventNumber << _format << _filename 
     << _unitchoice << _geneventPrecision;
}

void HIHepMCFile::persistentInput(PersistentIStream & is, int) {
  is >> _eventNumber >> _format >> _filename 
     >> _unitchoice >> _geneventPrecision;
}


ClassDescription<HIHepMCFile> HIHepMCFile::initHIHepMCFile;
// Definition of the static class description member.

void HIHepMCFile::Init() {

  static ClassDocumentation<HIHepMCFile> documentation
    ("This analysis handler will output the event record in HepMC format.");

  static Parameter<HIHepMCFile,long> interfacePrintEvent
    ("PrintEvent",
     "The number of events that should be printed.",
     &HIHepMCFile::_eventNumber, 1, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<HIHepMCFile,int> interfaceFormat
    ("Format",
     "Output format (1 = GenEvent, 2 = AsciiParticles, 5 = HepMC dump)",
     &HIHepMCFile::_format, 1, false, false);
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

  static Parameter<HIHepMCFile,string> interfaceFilename
    ("Filename", "Name of the output file",
     &HIHepMCFile::_filename, "");

  static Parameter<HIHepMCFile,unsigned int> interfacePrecision
    ("Precision",
     "Choice of output precision for the GenEvent format "
     " (as number of digits).",
     &HIHepMCFile::_geneventPrecision, 16, 6, 16,
     false, false, Interface::limited);
  
  static Switch<HIHepMCFile,int> interfaceUnits
    ("Units",
     "Unit choice for energy and length",
     &HIHepMCFile::_unitchoice, 0, false, false);
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
