// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLORivetAnalysis class.
//

#include <config.h>
#include "NLORivetAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Vectors/HepMCConverter.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"

using namespace ThePEG;

NLORivetAnalysis::NLORivetAnalysis() 
  :  _remnantId(82),_unitchoice(),
   debug(false), _rivet(), _nevent(0) {}

namespace {

// Special anonymous function for creating a genEvent.
HepMC::GenEvent * makeEvent(tEventPtr event, tSubProPtr sub, long no, long remnantId,
					  Energy eUnit, Length lUnit, 
					  CrossSection xsec, CrossSection xsecErr) {
  
  // generate beam particles
  const PPair& beam = event->incoming();
  HepMC::GenParticlePtr b1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum(),beam.first->id(),
					      1,eUnit);
  HepMC::GenParticlePtr b2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum(),beam.second->id(),
					      1,eUnit);

  // generate remnants
  HepMC::GenParticlePtr r1 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.first->momentum() - 
					      sub->incoming().first->momentum(),
					      remnantId,1,eUnit);
  HepMC::GenParticlePtr r2 =
    HepMCTraits<HepMC::GenEvent>::newParticle(beam.second->momentum() - 
					      sub->incoming().second->momentum(),
					      remnantId,1,eUnit);

  // generate outgoing particles
  vector<HepMC::GenParticlePtr> outgoing;
  for ( ParticleVector::const_iterator p = sub->outgoing().begin();
	p != sub->outgoing().end(); ++p ) {
    outgoing.push_back(HepMCTraits<HepMC::GenEvent>::newParticle((**p).momentum(),(**p).id(),
								 1,eUnit));
  }

  // generate one blob vertex
  HepMC::GenVertexPtr vertex = HepMCTraits<HepMC::GenEvent>::newVertex();

  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b1);
  HepMCTraits<HepMC::GenEvent>::addIncoming(*vertex,b2);

  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r1);
  HepMCTraits<HepMC::GenEvent>::addOutgoing(*vertex,r2);

  for ( vector<HepMC::GenParticlePtr>::const_iterator p = outgoing.begin();
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
}

void NLORivetAnalysis::analyze(ThePEG::tEventPtr event, long ieve, int loop, int state) {
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

  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // convert to hepmc

  HepMC::GenEvent * hepmc = 
    makeEvent(event,sub,_nevent,_remnantId,eUnit,lUnit,xsec,xsecErr);
 

  CurrentGenerator::Redirect stdout(cout);
  if(_rivet){
#if ThePEG_RIVET_VERSION == 1
    _rivet->analyze(const_cast<const HepMC::GenEvent&>(*hepmc));
#elif ThePEG_RIVET_VERSION > 1
    try {
      _rivet->analyze(const_cast<const HepMC::GenEvent&>(*hepmc));
    } catch (const YODA::Exception & e) {
      Throw<Exception>() << "Warning: Rivet/Yoda got the exception: "<< e.what()<<"\n"
                         << Exception::warning;
    }
#else
#error "Unknown ThePEG_RIVET_VERSION"
#endif
  }
  // delete hepmc event
  delete hepmc;
  
  if ( grp ) {

    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {

      hepmc = makeEvent(event,*s,_nevent,_remnantId,eUnit,lUnit,xsec,xsecErr);

      if ( _rivet ){
#if ThePEG_RIVET_VERSION == 1
        _rivet->analyze(const_cast<const HepMC::GenEvent&>(*hepmc));
#elif ThePEG_RIVET_VERSION > 1
        try {
          _rivet->analyze(const_cast<const HepMC::GenEvent&>(*hepmc));
        } catch (const YODA::Exception & e) {
          Throw<Exception>() << "Warning: Rivet/Yoda got the exception: "<< e.what()<<"\n"
                             << Exception::warning;
        }
#else
#error "Unknown ThePEG_RIVET_VERSION"
#endif
      }
      // delete hepmc event
      delete hepmc;
      

    }

  }

  ++_nevent;

}


ThePEG::IBPtr NLORivetAnalysis::clone() const {
  return new_ptr(*this);
}

ThePEG::IBPtr NLORivetAnalysis::fullclone() const {
  return new_ptr(*this);
}

void NLORivetAnalysis::persistentOutput(ThePEG::PersistentOStream & os) const {
  os << _analyses << filename << debug;
}

void NLORivetAnalysis::persistentInput(ThePEG::PersistentIStream & is, int) {
  is >> _analyses >> filename >> debug;
}

ThePEG::ClassDescription<NLORivetAnalysis> NLORivetAnalysis::initNLORivetAnalysis;
// Definition of the static class description member.

void NLORivetAnalysis::Init() {

  static ThePEG::ClassDocumentation<NLORivetAnalysis> documentation
    ("The NLORivetAnalysis class is a simple class to allow analyses"
     " from the Rivet library to be called from ThePEG");

  static ThePEG::ParVector<NLORivetAnalysis,string> interfaceAnalyses
    ("Analyses",
     "The names of the Rivet analyses to use",
     &NLORivetAnalysis::_analyses, -1, "", "","" "",
     false, false, ThePEG::Interface::nolimits);

  static Parameter<NLORivetAnalysis,long> interfaceRemnantId
    ("RemnantId",
     "Set the PDG id to be used for remnants.",
     &NLORivetAnalysis::_remnantId, 82, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<NLORivetAnalysis,string> interfaceFilename
    ("Filename",
#if ThePEG_RIVET_VERSION == 1
     "The name of the file where the AIDA histograms are put. If empty, "
     "the run name will be used instead. '.aida' will in any case be "
     "appended to the file name.",
#elif ThePEG_RIVET_VERSION > 1
     "The name of the file where the YODA histograms are put. If empty, "
     "the run name will be used instead. '.yoda' will in any case be "
     "appended to the file name.",
#else
#error "Unknown ThePEG_RIVET_VERSION"
#endif
     &NLORivetAnalysis::filename, "", true, false);


  static Switch<NLORivetAnalysis,bool> interfaceDebug
    ("Debug",
     "Enable debug information from Rivet",
     &NLORivetAnalysis::debug, false, true, false);
  static SwitchOption interfaceDebugNo
    (interfaceDebug,
     "No",
     "Disable debug information.",
     false);
  static SwitchOption interfaceDebugYes
    (interfaceDebug,
     "Yes",
     "Enable debug information from Rivet.",
     true);


  interfaceAnalyses.rank(10);

}

void NLORivetAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if( _nevent > 0 && _rivet ) {
    CurrentGenerator::Redirect stdout(cout);
#if ThePEG_RIVET_VERSION > 2
    _rivet->setCrossSection(make_pair(generator()->integratedXSec()/picobarn,
                                      generator()->integratedXSecErr()/picobarn));
#else
    _rivet->setCrossSection(generator()->integratedXSec()/picobarn);
#endif
    _rivet->finalize();

    string fname = filename;
#if ThePEG_RIVET_VERSION == 1
    if ( fname.empty() )
      fname = generator()->path() + "/" + generator()->runName() + ".aida";
#elif ThePEG_RIVET_VERSION > 1
    if ( fname.empty() )
      fname = generator()->path() + "/" + generator()->runName() + ".yoda";
#else
#error "Unknown ThePEG_RIVET_VERSION"
#endif
    _rivet->writeData(fname);
  }
  delete _rivet;
  _rivet = nullptr;
}

void NLORivetAnalysis::doinit() {
  AnalysisHandler::doinit();
  if(_analyses.empty()) 
    throw ThePEG::Exception() << "Must have at least one analysis loaded in "
			      << "NLORivetAnalysis::doinitrun()"
			      << ThePEG::Exception::runerror;

  // check that analysis list is available
  _rivet = new Rivet::AnalysisHandler; //(fname);
  _rivet->addAnalyses(_analyses);
  if ( _rivet->analysisNames().size() != _analyses.size() ) {
    throw ThePEG::Exception() 
      << "Rivet could not find all requested analyses.\n"
      << "Use 'rivet --list-analyses' to check availability.\n"
      << ThePEG::Exception::runerror;
  }
  delete _rivet;
  _rivet = 0;
}

void NLORivetAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // create NLORivet analysis handler
  CurrentGenerator::Redirect stdout(cout);
  _rivet = new Rivet::AnalysisHandler; //(fname);
  _rivet->addAnalyses(_analyses);
  // check that analysis list is still available
  if ( _rivet->analysisNames().size() != _analyses.size() ) {
    throw ThePEG::Exception() 
      << "Rivet could not find all requested analyses.\n"
      << "Use 'rivet --list-analyses' to check availability.\n"
      << ThePEG::Exception::runerror;
  }
  if ( debug )
    Rivet::Log::setLevel("Rivet",Rivet::Log::DEBUG);
}
