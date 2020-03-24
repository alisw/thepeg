// -*- C++ -*-
//
// HIHepMCFile.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HIHepMCFile class.
//

#include "HIHepMCFile.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Vectors/HepMCConverter.h"

using namespace ThePEG;

HIHepMCFile::HIHepMCFile() 
  : _eventNumber(1), _format(1), _filename(),
#ifdef HAVE_HEPMC_ROOTIO 
    _ttreename(),_tbranchname(),
#endif
    _unitchoice(), _geneventPrecision(16) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
HIHepMCFile::HIHepMCFile(const HIHepMCFile & x) 
  : AnalysisHandler(x), 
    _eventNumber(x._eventNumber), _format(x._format), 
    _filename(x._filename), 
#ifdef HAVE_HEPMC_ROOTIO 
    _ttreename(x._ttreename),_tbranchname(x._tbranchname),
#endif        
    _hepmcio(), _hepmcdump(), 
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

  switch ( _format ) {
#ifdef HAVE_HEPMC3
  default: {
    if ( _filename.empty() )
      _filename = generator()->filename() + ".hepmc";
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
      = new HepMC::WriterAscii(_filename.c_str());
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
  case 8:   {
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

void HIHepMCFile::dofinish() {
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
#ifdef HAVE_HEPMC3
    if (!_hepmcio->run_info()) 
    {
       _hepmcio->set_run_info(std::make_shared<HepMC::GenRunInfo>());
    std::vector<std::string>  w_names;
    w_names.push_back("Default");
    for ( map<string,double>::const_iterator w = event->optionalWeights().begin();
     w != event->optionalWeights().end(); ++w ) {
     w_names.push_back(w->first);
    }
    _hepmcio->run_info()->set_weight_names(w_names);  
    }
#endif
  HepMC::GenEvent * hepmc 
    = HepMCConverter<HepMC::GenEvent>::convert(*event, false,
					       eUnit, lUnit);
#ifdef HAVE_HEPMC3  
  std::shared_ptr<HepMC::HeavyIon> heavyion=std::make_shared<HepMC::HeavyIon>();
  heavyion->set(1,1,1,1,1,1);
#else
  HepMC::HeavyIon heavyion(1,1,1,1,1,1);
#endif  
  const LorentzPoint v1 = event->incoming().first->vertex();
  const LorentzPoint v2 = event->incoming().second->vertex();  
  double bpar = (v1 - v2).perp()/femtometer;

#ifdef HAVE_HEPMC3  
  heavyion->event_plane_angle=atan2((v1 - v2).y(),(v1 - v2).x());
  heavyion->impact_parameter=float(bpar);
#else
  heavyion.HepMC::HeavyIon::set_event_plane_angle(atan2((v1 - v2).y(),(v1 - v2).x()));
  heavyion.HepMC::HeavyIon::set_impact_parameter(float(bpar));  
#endif


#ifdef HAVE_HEPMC3  
  // Clear and blatant abuse of the Pdf info container!!
  HepMC::GenPdfInfoPtr pdfinfo=std::make_shared<HepMC::GenPdfInfo>();
  pdfinfo->set(1,1,event->optionalWeight("averageKappa"),event->optionalWeight("junctions"),event->optionalWeight("lambdaSum"),1,1);
#else
  HepMC::PdfInfo pdfinfo(1,1,event->optionalWeight("averageKappa"),event->optionalWeight("junctions"),event->optionalWeight("lambdaSum"),1,1);
#endif
  
#ifdef HAVE_HEPMC3  
  hepmc->set_heavy_ion(heavyion);
  hepmc->set_pdf_info(pdfinfo);
#else
  hepmc->set_heavy_ion(heavyion);
  hepmc->set_pdf_info(pdfinfo);
#endif

#ifdef HAVE_HEPMC3
    if (!_hepmcio->run_info())
    {
    _hepmcio->set_run_info(std::make_shared<HepMC::GenRunInfo>());
    std::vector<std::string>  w_names;
    w_names.push_back("Default");
    for ( map<string,double>::const_iterator w = event->optionalWeights().begin();
     w != event->optionalWeights().end(); ++w ) {
     w_names.push_back(w->first);
    }
    _hepmcio->run_info()->set_weight_names(w_names);  
    }
    hepmc->set_run_info(_hepmcio->run_info());
  if (_hepmcio)
    _hepmcio->write_event(*hepmc);
#else
  if (_hepmcio)
    _hepmcio->write_event(hepmc);
  else
    hepmc->print(_hepmcdump);
#endif
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
#ifdef HAVE_HEPMC3
#ifdef HAVE_HEPMC_ROOTIO    
     "Output format (1 = GenEvent,  6 = GenEventHepMC3, 7 = HEPEVT, 8 = GenEvent in ROOT, 9 = GenEvent in ROOT TTree  )",
#else
     "Output format (1 = GenEvent,  6 = GenEventHepMC3, 7 = HEPEVT",
#endif
#else
     "Output format (1 = GenEvent, 2 = AsciiParticles, 5 = HepMC dump",
#endif
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
#ifdef HAVE_HEPMC_ROOTIO 
  static Parameter<HIHepMCFile,string> interfaceTTreename
    ("TTreename", "Name of the TTree in output file",
     &HIHepMCFile::_ttreename, "hepmc3_tree");  
  static Parameter<HIHepMCFile,string> interfaceTBranchname
    ("TBranchname", "Name of the branch in output file",
     &HIHepMCFile::_tbranchname, "hepmc3_tree");
#endif
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
