// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOHepMC3File class.
//

#include "NLOHepMC3File.h"
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
#include "HepMC/IO/IO_GenEvent.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

NLOHepMC3File::NLOHepMC3File() 
  : _remnantId(82), _filename(), _unitchoice(),
    _geneventPrecision(16), _eventNumber(1) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
NLOHepMC3File::NLOHepMC3File(const NLOHepMC3File & x) 
  : AnalysisHandler(x), 
    _remnantId(x._remnantId), 
    _filename(x._filename), _hepmcio(),
    _unitchoice(x._unitchoice), 
    _geneventPrecision(x._geneventPrecision),
    _eventNumber(x._eventNumber) {}

HepMC::GenEvent * NLOHepMC3File::makeEvent(tEventPtr event, tSubProPtr sub, long no,
					  Energy eUnit, Length lUnit, 
					  CrossSection xsec, CrossSection xsecErr) const {
  
  // generate beam particles
  const PPair& beam = event->incoming();
  HepMC::GenParticle * b1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum(),beam.first->id(),
					      1,eUnit);
  HepMC::GenParticle * b2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum(),beam.second->id(),
					      1,eUnit);

  // generate remnants
  HepMC::GenParticle * r1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum() - 
					      sub->incoming().first->momentum(),
					      _remnantId,1,eUnit);
  HepMC::GenParticle * r2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum() - 
					      sub->incoming().second->momentum(),
					      _remnantId,1,eUnit);

  // generate outgoing particles
  vector<HepMC::GenParticle*> outgoing;
  for ( ParticleVector::const_iterator p = sub->outgoing().begin();
	p != sub->outgoing().end(); ++p ) {
    outgoing.push_back(HepMCTraits<HepMC::GenEvent>::newParticle((**p).momentum(),(**p).id(),
								 1,eUnit));
  }

  // generate one blob vertex
  HepMC::GenVertex * vertex = HepMCTraits<HepMC::GenEvent>::newVertex();

  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b1);
  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b2);

  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r1);
  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r2);

  for ( vector<HepMC::GenParticle*>::const_iterator p = outgoing.begin();
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

void NLOHepMC3File::analyze(tEventPtr event, long, int, int) {

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

  HepMC::GenEvent * hepmc = 
    makeEvent(event,sub,_eventNumber,eUnit,lUnit,xsec,xsecErr);
  if (_hepmcio)
    _hepmcio->write_event(*hepmc);
  delete hepmc;

  if ( grp ) {

    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {

      hepmc = makeEvent(event,*s,_eventNumber,eUnit,lUnit,xsec,xsecErr);
      if (_hepmcio)
	_hepmcio->write_event(*hepmc);
      delete hepmc;

    }

  }

  ++_eventNumber;

}

IBPtr NLOHepMC3File::clone() const {
  return new_ptr(*this);
}

IBPtr NLOHepMC3File::fullclone() const {
  return new_ptr(*this);
}

void NLOHepMC3File::doinitrun() {
  AnalysisHandler::doinitrun();

  // set default filename unless user-specified name exists
  if ( _filename.empty() )
    _filename = generator()->filename() + ".hepmc";

  HepMC::IO_GenEvent * tmpio 
    = new HepMC::IO_GenEvent(_filename.c_str(), ios::out);
  tmpio->set_precision(_geneventPrecision);
  _hepmcio = tmpio;

}

void NLOHepMC3File::dofinish() {
  if (_hepmcio) {
    delete _hepmcio;
    _hepmcio = 0;
  }
  AnalysisHandler::dofinish();
  cout << "\nNLOHepMC3File: generated HepMC output.\n";
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOHepMC3File::persistentOutput(PersistentOStream & os) const {
  os << _remnantId << _filename 
     << _unitchoice << _geneventPrecision << _eventNumber;
}

void NLOHepMC3File::persistentInput(PersistentIStream & is, int) {
  is >> _remnantId >> _filename 
     >> _unitchoice >> _geneventPrecision >> _eventNumber;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOHepMC3File,AnalysisHandler>
  describeHerwigNLOHepMC3File("ThePEG::NLOHepMCFile", "HepMCAnalysis.so");

void NLOHepMC3File::Init() {

  static ClassDocumentation<NLOHepMC3File> documentation
    ("Write hard sub processes or sub process groups to HepMC.");


  static Parameter<NLOHepMC3File,long> interfaceRemnantId
    ("RemnantId",
     "Set the PDG id to be used for remnants.",
     &NLOHepMC3File::_remnantId, 82, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<NLOHepMC3File,string> interfaceFilename
    ("Filename", "Name of the output file",
     &NLOHepMC3File::_filename, "");

  static Parameter<NLOHepMC3File,unsigned int> interfacePrecision
    ("Precision",
     "Choice of output precision for the GenEvent format "
     " (as number of digits).",
     &NLOHepMC3File::_geneventPrecision, 16, 6, 16,
     false, false, Interface::limited);
  
  static Switch<NLOHepMC3File,int> interfaceUnits
    ("Units",
     "Unit choice for energy and length",
     &NLOHepMC3File::_unitchoice, 0, false, false);
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

