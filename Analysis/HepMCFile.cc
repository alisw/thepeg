// -*- C++ -*-
//
// HepMCFile.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HepMCFile class.
//

#include "HepMCFile.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Vectors/HepMCConverter.h"
using namespace ThePEG;

HepMCFile::HepMCFile() 
  : _eventNumber(1), _format(1), _filename(),
#ifdef HAVE_HEPMC_ROOTIO 
   _ttreename(),_tbranchname(),
#endif
    _unitchoice(), _geneventPrecision(16), _addHI(0) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
HepMCFile::HepMCFile(const HepMCFile & x) 
  : AnalysisHandler(x), 
    _eventNumber(x._eventNumber), _format(x._format), 
    _filename(x._filename),
#ifdef HAVE_HEPMC_ROOTIO 
    _ttreename(x._ttreename),_tbranchname(x._tbranchname),
#endif
    _hepmcio(), _hepmcdump(), _unitchoice(x._unitchoice), 
    _geneventPrecision(x._geneventPrecision) {}

IBPtr HepMCFile::clone() const {
  return new_ptr(*this);
}

IBPtr HepMCFile::fullclone() const {
  return new_ptr(*this);
}

void HepMCFile::doinitrun() {
  AnalysisHandler::doinitrun();

   if ( _filename.empty() )
      _filename = generator()->filename() + ".hepmc";


  switch ( _format ) {
#ifdef HAVE_HEPMC3
  default: {
    HepMC::WriterAsciiHepMC2 * tmpio 
      = new HepMC::WriterAsciiHepMC2(_filename.c_str());
    tmpio->set_precision(_geneventPrecision);
    _hepmcio = tmpio;
  }
    break;
  case 6: {
    if ( _filename.empty() )
      _filename = generator()->filename() + ".hepmc";
    HepMC::WriterAscii * tmpio 
      = new HepMC::WriterAscii(_filename.c_str(),NULL);
    tmpio->set_precision(_geneventPrecision);
    _hepmcio = tmpio;
  }
    break;
  case 7: {
    if ( _filename.empty() )
      _filename = generator()->filename() + ".hepevt";
    HepMC::WriterHEPEVT * tmpio 
      = new  HepMC::WriterHEPEVT(_filename.c_str()); 
    _hepmcio = tmpio;
  }
    break;
#ifdef HAVE_HEPMC_ROOTIO
  case 8: {
    if ( _filename.empty() )
      _filename = generator()->filename() + ".root";
    HepMC::WriterRoot * tmpio 
      = new HepMC::WriterRoot(_filename.c_str());
    _hepmcio = tmpio;  
  }
    break;
  case 9:   {
    if ( _filename.empty() )
      _filename = generator()->filename() + ".root";
    HepMC::WriterRootTree * tmpio 
      = new HepMC::WriterRootTree(_filename.c_str());
    _hepmcio = tmpio;  
  }
    break;
#endif
#else
  default: {
    HepMC::IO_GenEvent * tmpio 
      = new HepMC::IO_GenEvent(_filename.c_str(), ios::out);
    tmpio->precision(_geneventPrecision);
    _hepmcio = tmpio;
  }
    break;
  case 2: 
    _hepmcio = new HepMC::IO_AsciiParticles(_filename.c_str(), ios::out); 
    break;
  case 5: 
    _hepmcio = 0; 
    _hepmcdump.open(_filename.c_str()); 
    break;
#endif
  }
}
void HepMCFile::dofinish() {
#ifdef HAVE_HEPMC3
  _hepmcio->close();
  delete _hepmcio;
#else
  if (_hepmcio) {
    delete _hepmcio;
    _hepmcio = 0;
  }
  else
    _hepmcdump.close();
#endif
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
#ifdef HAVE_HEPMC3
    _hepmcio->set_run_info(std::make_shared<HepMC::GenRunInfo>());
    std::vector<std::string>  w_names;
    w_names.push_back("Default");
    for ( map<string,double>::const_iterator w = event->optionalWeights().begin();
     w != event->optionalWeights().end(); ++w ) {
     w_names.push_back(w->first);
    }
    _hepmcio->run_info()->set_weight_names(w_names);  
#endif

  HepMC::GenEvent * hepmc 
    = HepMCConverter<HepMC::GenEvent>::convert(*event, false,
					       eUnit, lUnit);

  const LorentzPoint v1 = event->incoming().first->vertex();
  const LorentzPoint v2 = event->incoming().second->vertex();  
  if (  _addHI > 0 || ( _addHI == 0 && v1.perp() >= ZERO && v2.perp() >= ZERO ) ) {
  double bpar = (v1 - v2).perp()/femtometer;

#ifdef HAVE_HEPMC3  
  std::shared_ptr<HepMC::HeavyIon> heavyion=std::make_shared<HepMC::HeavyIon>();
  heavyion->set(1,1,1,1,1,1);
  heavyion->event_plane_angle=atan2((v1 - v2).y(),(v1 - v2).x());
  heavyion->impact_parameter=float(bpar);
  hepmc->set_heavy_ion(heavyion);
#else
  HepMC::HeavyIon heavyion(1,1,1,1,1,1);
  heavyion.HepMC::HeavyIon::set_event_plane_angle(atan2((v1 - v2).y(),(v1 - v2).x()));
  heavyion.HepMC::HeavyIon::set_impact_parameter(float(bpar));  
#endif  

  hepmc->set_heavy_ion(heavyion);

  }

#ifdef HAVE_HEPMC3
  hepmc->set_run_info( _hepmcio->run_info()); 
  _hepmcio->write_event(*hepmc);
#else
  if (_hepmcio)
    _hepmcio->write_event(hepmc);
  else
    hepmc->print(_hepmcdump);
#endif

  delete hepmc;

}

void HepMCFile::persistentOutput(PersistentOStream & os) const {
  os << _eventNumber << _format << _filename 
     << _unitchoice << _geneventPrecision << _addHI;
}

void HepMCFile::persistentInput(PersistentIStream & is, int) {
  is >> _eventNumber >> _format >> _filename 
     >> _unitchoice >> _geneventPrecision >> _addHI;
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
#ifdef HAVE_HEPMC3
#ifdef HAVE_HEPMC_ROOTIO    
     "Output format (1 = GenEvent,  6 = GenEventHepMC3, 7 = HEPEVT, 8 = GenEvent in ROOT, 9 = GenEvent in ROOT TTree  )",
#else
     "Output format (1 = GenEvent,  6 = GenEventHepMC3, 7 = HEPEVT",
#endif
#else
     "Output format (1 = GenEvent, 2 = AsciiParticles, 5 = HepMC dump",
#endif
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
#ifdef HAVE_HEPMC3
  static SwitchOption interfaceFormatGenEventHepMC3
    (interfaceFormat,
     "GenEventHepMC3",
     "GenEvent in HepMC3",
     6);
  static SwitchOption interfaceFormatHEPEVT
    (interfaceFormat,
     "HEPEVT",
     "HEPEVT",
     7);
#ifdef HAVE_HEPMC_ROOTIO 
  static SwitchOption interfaceFormatGenEventROOT
    (interfaceFormat,
     "GenEventROOT",
     "GenEvent in ROOT",
     8);
  static SwitchOption interfaceFormatGenEventROOTTree
    (interfaceFormat,
     "GenEventROOTTree",
     "GenEvent in ROOT TTree",
     9);     
#endif
#endif



  static Parameter<HepMCFile,string> interfaceFilename
    ("Filename", "Name of the output file",
     &HepMCFile::_filename, "");
#ifdef HAVE_HEPMC_ROOTIO 
  static Parameter<HepMCFile,string> interfaceTTreename
    ("TTreename", "Name of the TTree in output file",
     &HepMCFile::_ttreename, "hepmc3_tree");  
  static Parameter<HepMCFile,string> interfaceTBranchname
    ("TBranchname", "Name of the branch in output file",
     &HepMCFile::_tbranchname, "hepmc3_tree");
#endif
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

  static Switch<HepMCFile,int> interfaceAddHI
    ("AddHI",
     "Options for adding heavy ion info to GenEvent.",
     &HepMCFile::_addHI, 0, false, false);
  static SwitchOption interfaceAddHIMaybe
    (interfaceAddHI,
     "Maybe",
     "Add Heavy Ion info if both incoming particles impact parameter is not exactly zero.",
     0);
  static SwitchOption interfaceAddHINever
    (interfaceAddHI,
     "Never",
     "Never add Heavy Ion info.",
     -1);
  static SwitchOption interfaceAddHIAlways
    (interfaceAddHI,
     "Always",
     "Always add Heavy Ion info.",
     1);

}
