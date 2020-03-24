// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOHepMCFile class.
//

#include "NLOHepMCFile.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Vectors/HepMCConverter.h"
using namespace ThePEG;

NLOHepMCFile::NLOHepMCFile() 
  : _remnantId(82), _format(1), _filename(), 
#ifdef HAVE_HEPMC_ROOTIO 
    _ttreename(),_tbranchname(),
#endif  
    _unitchoice(),_geneventPrecision(16), _eventNumber(1) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
NLOHepMCFile::NLOHepMCFile(const NLOHepMCFile & x) 
  : AnalysisHandler(x), 
    _remnantId(x._remnantId), _format(x._format), 
    _filename(x._filename),
#ifdef HAVE_HEPMC_ROOTIO 
    _ttreename(x._ttreename),_tbranchname(x._tbranchname),
#endif
    _hepmcio(), _hepmcdump(), _unitchoice(x._unitchoice), 
    _geneventPrecision(x._geneventPrecision),
    
    _eventNumber(x._eventNumber) {}

HepMC::GenEvent * NLOHepMCFile::makeEvent(tEventPtr event, tSubProPtr sub, long no,
					  Energy eUnit, Length lUnit, 
					  CrossSection xsec, CrossSection xsecErr) const {

  typedef typename HepMCTraits<HepMC::GenEvent>::ParticlePtrT ParticlePtrT;
  typedef typename HepMCTraits<HepMC::GenEvent>::VertexPtrT VertexPtrT;
  
  // generate beam particles
  const PPair& beam = event->incoming();
  ParticlePtrT b1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum(),beam.first->id(),
					      1,eUnit);
  ParticlePtrT b2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum(),beam.second->id(),
					      1,eUnit);

  // generate remnants
  ParticlePtrT r1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum() - 
					      sub->incoming().first->momentum(),
					      _remnantId,1,eUnit);
  ParticlePtrT r2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum() - 
					      sub->incoming().second->momentum(),
					      _remnantId,1,eUnit);

  // generate outgoing particles
  vector<ParticlePtrT> outgoing;
  for ( ParticleVector::const_iterator p = sub->outgoing().begin();
	p != sub->outgoing().end(); ++p ) {
    outgoing.push_back(HepMCTraits<HepMC::GenEvent>::newParticle((**p).momentum(),(**p).id(),
								 1,eUnit));
  }

  // generate one blob vertex
  VertexPtrT vertex = HepMCTraits<HepMC::GenEvent>::newVertex();

  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b1);
  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b2);

  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r1);
  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r2);

  for ( vector<ParticlePtrT>::const_iterator p = outgoing.begin();
	p != outgoing.end(); ++p )
    HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,*p);

  HepMC::GenEvent * ev = 
    HepMCTraits<HepMC::GenEvent>::newEvent(no,event->weight()*sub->groupWeight(),
					   event->optionalWeights());
  HepMCTraits<HepMC::GenEvent>::setUnits(*ev,eUnit,lUnit);
  HepMCTraits<HepMC::GenEvent>::setBeamParticles(*ev,b1,b2);

  HepMCTraits<HepMC::GenEvent>::addVertex(*ev,vertex);

  HepMCTraits<HepMC::GenEvent>::setCrossSection(*ev,xsec/picobarn,
						xsecErr/picobarn);

  return ev;

}

void NLOHepMCFile::analyze(tEventPtr event, long, int, int) {

  Energy eUnit;
  Length lUnit;
  switch (_unitchoice) {
  default: eUnit = GeV; lUnit = millimeter; break;
  case 1:  eUnit = MeV; lUnit = millimeter; break;
  case 2:  eUnit = GeV; lUnit = centimeter; break;
  case 3:  eUnit = MeV; lUnit = centimeter; break;
  }

  tcEHPtr eh = dynamic_ptr_cast<tcEHPtr>(event->primaryCollision()->handler());
  assert(eh);

  CrossSection xsec = eh->integratedXSec();
  CrossSection xsecErr = eh->integratedXSecErr();

  tSubProPtr sub = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr grp = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);
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
  HepMC::GenEvent * hepmc = 
    makeEvent(event,sub,_eventNumber,eUnit,lUnit,xsec,xsecErr);
#ifdef HAVE_HEPMC3
    _hepmcio->write_event(*hepmc);
#else
  if (_hepmcio)
    _hepmcio->write_event(hepmc);
  else
    hepmc->print(_hepmcdump);
#endif
  delete hepmc;

  if ( grp ) {

    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {

      hepmc = makeEvent(event,*s,_eventNumber,eUnit,lUnit,xsec,xsecErr);
#ifdef HAVE_HEPMC3
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

  }

  ++_eventNumber;

}

IBPtr NLOHepMCFile::clone() const {
  return new_ptr(*this);
}

IBPtr NLOHepMCFile::fullclone() const {
  return new_ptr(*this);
}

void NLOHepMCFile::doinitrun() {
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
      = new HepMC::WriterRootTree(_filename.c_str(),_ttreename.c_str(),_tbranchname.c_str());
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

void NLOHepMCFile::dofinish() {
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
   cout << "\nNLOHepMCFile: generated HepMC output.\n";
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOHepMCFile::persistentOutput(PersistentOStream & os) const {
  os << _remnantId << _format << _filename 
     << _unitchoice << _geneventPrecision << _eventNumber;
}

void NLOHepMCFile::persistentInput(PersistentIStream & is, int) {
  is >> _remnantId >> _format >> _filename 
     >> _unitchoice >> _geneventPrecision >> _eventNumber;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOHepMCFile,AnalysisHandler>
  describeHerwigNLOHepMCFile("ThePEG::NLOHepMCFile", "HepMCAnalysis.so");

void NLOHepMCFile::Init() {

  static ClassDocumentation<NLOHepMCFile> documentation
    ("Write hard sub processes or sub process groups to HepMC.");


  static Parameter<NLOHepMCFile,long> interfaceRemnantId
    ("RemnantId",
     "Set the PDG id to be used for remnants.",
     &NLOHepMCFile::_remnantId, 82, 0, 0,
     false, false, Interface::nolimits);

  static Switch<NLOHepMCFile,int> interfaceFormat
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
     &NLOHepMCFile::_format, 1, false, false);
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

  static Parameter<NLOHepMCFile,string> interfaceFilename
    ("Filename", "Name of the output file",
     &NLOHepMCFile::_filename, "");
#ifdef HAVE_HEPMC_ROOTIO 
  static Parameter<NLOHepMCFile,string> interfaceTTreename
    ("TTreename", "Name of the TTree in output file",
     &NLOHepMCFile::_ttreename, "hepmc3_tree");  
  static Parameter<NLOHepMCFile,string> interfaceTBranchname
    ("TBranchname", "Name of the branch in output file",
     &NLOHepMCFile::_tbranchname, "hepmc3_tree");
#endif
  static Parameter<NLOHepMCFile,unsigned int> interfacePrecision
    ("Precision",
     "Choice of output precision for the GenEvent format "
     " (as number of digits).",
     &NLOHepMCFile::_geneventPrecision, 16, 6, 16,
     false, false, Interface::limited);
  
  static Switch<NLOHepMCFile,int> interfaceUnits
    ("Units",
     "Unit choice for energy and length",
     &NLOHepMCFile::_unitchoice, 0, false, false);
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

