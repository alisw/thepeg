// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerHandler class.
//

#include "ShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PDF.h"
#include "ThePEG/PDF/PDFBase.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Pythia7/Shower/Beam.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/Throw.h"

#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "ShowerHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Pythia7;

ShowerHandler::~ShowerHandler() {
  if ( theShowerModel ) delete theShowerModel;
}

Shower::Shower * ShowerHandler::getModel() {
  if ( !theShowerModel ) theShowerModel = new Shower::Shower;
  theShowerModel->ISR = 0;
  theShowerModel->theSpaceShower = 0;
  theShowerModel->FSR = 0;
  theShowerModel->theTimeShower = 0;
  theShowerModel->FSRONISR = 0;
  
  if ( timeShower() ) {
    theShowerModel->theTimeShower = timeShower()->getModel();
    theShowerModel->FSR = 1;
  }
  if ( spaceShower() ) {
    theShowerModel->theSpaceShower = spaceShower()->getModel();
    theShowerModel->ISR = 1;
  }
  if ( addFSROnISR  && timeShower() && spaceShower() )
    theShowerModel->FSRONISR = 1;
  return theShowerModel;
}

void ShowerHandler::
handle(EventHandler & ch, const tPVector & tagged, const Hint & hint) {

  tPVector final;

  if ( hint.tagged() && !tagged.empty() ) {

    for ( int i = 0, N = tagged.size(); i < N; ++i )
      if ( canShower(tagged[i]) ) final.push_back(tagged[i]);

    if ( !final.empty() ) cascade(final);

    return;

  }

  tCollPtr coll = StepHandler::eventHandler()->currentCollision();
  for ( int i = 0, N = coll->subProcesses().size(); i < N; ++i ) {
    if ( coll->subProcesses()[i]->decayed() ) continue;
    for ( int itry = 0; itry < maxTries; ++itry ) {
      try {
	cascade(coll->subProcesses()[i]);
	break;
      } catch (Veto) {}
      Throw<InfiniteLoopException>()
	<< "Shower handler '" << name() << "' giving up after trying "
	<< maxTries << " times to showe the sub-process."
	<< Exception::eventerror;
    }
  }

  for ( ParticleSet::iterator it = currentStep()->particles().begin();
	it != currentStep()->particles().end(); ++it )
    if ( canShower(*it) ) final.push_back(*it);

  cascade(final);

}

void ShowerHandler::cascade(const tPVector & final) {}

void ShowerHandler::cascade(tSubProPtr sub) {
  sub->decayed(true);
  colourIndex.clear();
  colourIndex(0, tColinePtr());
  particleIndex.clear();
  PartonExtractor & pex = *StepHandler::eventHandler()->lastExtractor();

  // First add parents if the incoming partons should be showered.
  event.zero();
  tPPair incoming = sub->incoming();
  tPPair beams;
  pair<Shower::PDF*,Shower::PDF*> pdfs;

  if ( canShower(incoming.first, true) ) {
    tPBIPtr partonBin = pex.partonBinInstance(incoming.first);
    if ( partonBin && partonBin->pdf() ) {
      beams.first = partonBin->particle();
      pdfs.first =
	new Shower::ThePEGPDF(partonBin->pdf(), beams.first->dataPtr());
    }
  }
  if ( canShower(incoming.second, true) ) {
    tPBIPtr partonBin = pex.partonBinInstance(incoming.second);
    if ( partonBin && partonBin->pdf() ) {
      beams.second = partonBin->particle();
      pdfs.second =
	new Shower::ThePEGPDF(partonBin->pdf(), beams.second->dataPtr());
    }
  }

  // First boost to cms.
  LorentzRotation r;
  if ( pdfs.first && pdfs.second ) {
    r = Utilities::boostToCM(beams);
    incoming.first->transform(r);
    incoming.second->transform(r);
 }
  else if ( pdfs.first && !pdfs.second ) {
    r = Utilities::boostToCM(make_pair(incoming.second, beams.first));
    incoming.first->transform(r);
  }    
  else if ( !pdfs.first && pdfs.second ) {
    r = Utilities::boostToCM(make_pair(incoming.first, beams.second));
    incoming.second->transform(r);
  }
  else {
    r = Utilities::boostToCM(make_pair(incoming.first, incoming.second));
  }

  Utilities::transform(sub->intermediates(), r);
  for ( int i = 0, N = sub->outgoing().size(); i < N; ++i )
    sub->outgoing()[i]->deepTransform(r);

  // Now add beam and incoming particles to the internal event.
  Shower::BeamParticle beam1;
  Shower::BeamParticle beam2;
  long mother1 = -1;
  long mother2 = -1;
  if ( pdfs.first && pdfs.second ) {
    addParticle(beams.first, -9, -1, -1);
    addParticle(beams.second, -9, -1, -1);
    addParticle(incoming.first, -1, 0, -1);
    addParticle(incoming.second, -1, 1, -1);
    beam1 = event[0];
    beam2 = event[1];
    beam1.pdf(pdfs.first);
    beam2.pdf(pdfs.second);
    mother1 = 2;
    mother2 = 3;
  }
  else if ( pdfs.first && !pdfs.second ) {
    addParticle(incoming.second, -1, -1, -1);
    addParticle(beams.first, -9, -1, -1);
    addParticle(incoming.first, -1, 1, -1);
    beam2 = event[1];
    beam2.pdf(pdfs.first);
    mother1 = 0;
    mother2 = 2;
  }
  else if ( !pdfs.first && pdfs.second ) {
    addParticle(incoming.first, -1, -1, -1);
    addParticle(beams.second, -9, -1, -1);
    addParticle(incoming.second, -1, 1, -1);
    beam2 = event[1];
    beam2.pdf(pdfs.first);
    mother1 = 0;
    mother2 = 2;
  }
  else {
    addParticle(incoming.first, -1, -1, -1);
    addParticle(incoming.second, -1, -1, -1);
    mother1 = 0;
    mother2 = 1;
  }

  // First add outgoing partons which do not come from resonances.
  for ( int i = 0, N = sub->outgoing().size(); i < N; ++i ) {
    if ( sub->outgoing()[i]->parents().size() != 1 ||
	 sub->outgoing()[i]->parents()[0] == incoming.first ||
	 sub->outgoing()[i]->parents()[0] == incoming.second ||
	 !isResonance(sub->outgoing()[i]->parents()[0]) ) {
      addFinalParticles(sub->outgoing()[i], mother1, mother2);
    }
  }

  // Then add outgoing partons which do come from resonances.
  for ( int i = 0, N = sub->outgoing().size(); i < N; ++i ) {
    if ( sub->outgoing()[i]->parents().size() == 1 &&
	 sub->outgoing()[i]->parents()[0] != incoming.first &&
	 sub->outgoing()[i]->parents()[0] != incoming.second &&
	 isResonance(sub->outgoing()[i]->parents()[0]) ) {
      long ri = addParticle(sub->outgoing()[i]->parents()[0],
			    2, mother1, mother2);
      addFinalParticles(sub->outgoing()[i], ri, -1);
    }
  }

  long outSize = event.size();

  //  CurrentGenerator::log() << event << endl;

  getModel()->shower(event, beam1, beam2, true);

  //  CurrentGenerator::log() << event;

  pair<long,long> inc1(-1, -1);
  pair<long,long> inc2(-1, -1);

  tPPair newincs = incoming;

  if ( pdfs.first && pdfs.second ) {
    findChain(inc1, 0);
    findChain(inc2, 1);
  }
  else if ( pdfs.first && !pdfs.second ) {
    findChain(inc1, 1);
  }
  else if ( !pdfs.first && pdfs.second ) {
    findChain(inc2, 1);
  }

  if ( pdfs.first ) {
    if ( event[inc1.second].status() == -2 )
      newincs.first = createParticle(inc1.first, true);
    else
      newincs.first->set5Momentum(momentum(inc1.first));
  }
  if ( pdfs.second ) {
    if ( event[inc2.second].status() == -2 )
      newincs.second = createParticle(inc2.first, true);
    else
      newincs.second->set5Momentum(momentum(inc2.first));
  }

  Lorentz5Momentum p1 = newincs.first->momentum();
  Lorentz5Momentum p2 = newincs.second->momentum();
  PBIPair newbins = pex.newRemnants(incoming, newincs, newStep());
  LorentzRotation tot = pex.boostRemnants(newbins, p1, p2,
					  pdfs.first, pdfs.second);

  for ( int i = outSize, N = event.size(); i < N; ++i ) {
    if ( particleIndex.included(i) ) continue;
    if ( !history() ) {
      if ( event[i].prev() >= 0 && event[i].prev() < outSize ) {
	if ( !particleIndex.included(event[i].prev()) )
	  event[i].prev(-1);
	else if ( event[i].status() != 1 ) 
	  particleIndex(i, getParticle(event[i].prev()));
      }
      if ( event[i].status() != 1 ) continue;
    }
    tPPtr p = createParticle(i);
    p->transform(tot);
    tcPPair pp = findParent(i);
    if ( pp.second ) {
      tcPVector parents;
      parents.push_back(pp.first);
      parents.push_back(pp.second);
      newStep()->addDecayProduct(parents.begin(), parents.end(), p);
    } else
      newStep()->addDecayProduct(pp.first, p, false);

    if ( i == inc1.second ) newStep()->setCopy(incoming.first, p);
    if ( i == inc2.second ) newStep()->setCopy(incoming.second, p);
    if ( event[i].prev() >= 0 && particleIndex.included(event[i].prev()) )
      newStep()->setCopy(getParticle(event[i].prev()), p);
  }

  if ( newincs.first != incoming.first )
    newincs.first->addChild(incoming.first);
  if ( newincs.second != incoming.second )
    newincs.second->addChild(incoming.second);

  if ( pdfs.first ) delete pdfs.first;
  if ( pdfs.second ) delete pdfs.second;

}

void ShowerHandler::
findChain(pair<long,long> & inc, long init) const {
  for ( int i = event.size() - 1; i > init; --i )
    if ( event[i].mother1() == init && event[i].mother2() < 0 &&
	 event[i].status() < 0 ) {
      if ( inc.first < 0 ) inc.first = i;
      else inc.second = i;
      findChain(inc, i);
      return;
    }
}

tPPtr ShowerHandler::getParticle(long i) const {
  return particleIndex(i);
}

tPPair ShowerHandler::findParent(long i) const {
  tPPair pp;
  if ( i < 0 ) return pp;
  pp = make_pair(particleIndex(event[i].mother1()),
		 particleIndex(event[i].mother2()));
  if ( pp.first == pp.second ) pp.second = tPPtr();
  return pp.first? pp: findParent(event[i].mother1());
}

Lorentz5Momentum ShowerHandler::momentum(long i) const {
  return Lorentz5Momentum(event[i].px()*GeV, event[i].py()*GeV,
			  event[i].pz()*GeV, event[i].e()*GeV,
			  event[i].m()*GeV);
}

tPPtr ShowerHandler::copyParticle(long i, tPPtr p) {
  p->set5Momentum(momentum(i));
  p->scale(sqr(event[i].scale()*GeV));
  tColinePtr c = colourIndex(event[i].col());
  if ( c ) c->addColoured(p);
  c = colourIndex(event[i].anticol());
  if ( c ) c->addAntiColoured(p);
  particleIndex(i, p);
  return p;
}

tPPtr ShowerHandler::createParticle(long i, bool inc) {
  tPPtr p = copyParticle(i, getParticleData(event[i].id())->produceParticle());
  if ( !inc && event[i].status() != 1 ) return p;
  if ( inc && !spaceShower() ) return p;
  if ( !inc && !timeShower() ) return p;
  if ( LeptonMatcher::Check(p->id()) )
    p->scale(sqr(inc? spaceShower()->Q0ChgL(): timeShower()->Q0ChgL()));
  if ( StandardQCDPartonMatcher::Check(p->id()) ||
       abs(p->id()) == ParticleID::t )
    p->scale(sqr(inc? spaceShower()->Q0(): timeShower()->Q0()));
  return p;
}

bool ShowerHandler::isResonance(tcPPtr r) const {
  if ( r->momentum().m2() <= 0.0*GeV2 ) return false;
  if ( r->children().size() < 2 ) return false;
  for ( int i = 0, N = r->children().size(); i < N; ++i ) {
    if ( r->children()[i]->parents().size() != 1 ) return false;
    if ( r->children()[i]->parents()[0] != r ) return false;
  }
  return true;
}

long ShowerHandler::addParticle(tPPtr p, long status) {
  long mother1 = -1;
  if ( p->parents().size() > 0 && particleIndex.included(p->parents()[0]) )
    mother1 = particleIndex(p->parents()[0]);
  long mother2 = -1;
  if ( p->parents().size() > 1 && particleIndex.included(p->parents()[1]) )
    mother2 = particleIndex(p->parents()[1]);
  return addParticle(p, status, mother1, mother2);
}

long ShowerHandler::
addParticle(tPPtr p, long status,
	    long mother1, long mother2) {
  if ( particleIndex.included(p) ) return particleIndex(p);
  long idx = event.size();
  particleIndex(idx, p);
  event.append(Shower::Particle(p->id(), status, mother1, mother2,
				colourIndex(p->colourLine()),
				colourIndex(p->antiColourLine()),
				p->momentum().x()/GeV, p->momentum().y()/GeV,
				p->momentum().z()/GeV, p->momentum().e()/GeV,
				p->momentum().mass()/GeV,
				sqrt(max(p->scale(), 0.0*GeV2))/GeV));
  return idx;
}

void ShowerHandler::
addFinalParticles(tPPtr p, long mother1, long mother2) {
  p = p->final();
  long status = 1;
  if ( p->decayed() ) status = 2;
  mother1 = addParticle(p, status, mother1, mother2);
  if ( status == 1 ) return;
  for ( int i = 0, N = p->children().size(); i < N; ++i )
    addFinalParticles(p->children()[i], mother1, -1);
}
  

bool ShowerHandler::canShower(tcPPtr p, bool inc) const {
  if ( inc && !spaceShower() ) return false;
  if ( !inc && !timeShower() ) return false;
  if ( LeptonMatcher::Check(p->id()) )
    return p->scale() > sqr(inc? spaceShower()->Q0ChgL():
			    timeShower()->Q0ChgL());
  if ( StandardQCDPartonMatcher::Check(p->id()) ||
       abs(p->id()) == ParticleID::t )
    return p->scale() > sqr(inc? spaceShower()->Q0(): timeShower()->Q0());

  return false;
}



void ShowerHandler::cascade() {}

void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << theTimeShower << theSpaceShower << addFSROnISR << theFullHistory
     << maxTries;
}

void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> theTimeShower >> theSpaceShower >> addFSROnISR >> theFullHistory
     >> maxTries;
}

ClassDescription<ShowerHandler> ShowerHandler::initShowerHandler;
// Definition of the static class description member.

void ShowerHandler::Init() {

  static ClassDocumentation<ShowerHandler> documentation
    ("The Pythia7::ShowerHandler class administers parton showers for "
     "external partons in a hard ThePEG::SubProcess or in a "
     "decay of a heavy resonance.",
     "Parton cascades performed according to the Pythia scheme for initial- "
     "\\cite{Sjo85,Ben85,Miu99} and final-state "
     "\\cite{Ben87,Nor01} parton showers.",
     "\\bibitem{Ben87} "
     "M.~Bengtsson and T.~Sj\\\"ostrand, Phys.~Lett.~{\\bf B185} (1987) 435;\n"
     "Nucl. Phys. B289 (1987) 810\n"
     "\\bibitem{Nor01} "
     "E.~Norrbin and T.~Sj\\\"ostrand, Nucl.~Phys.{\\bf B603} (2001) 297\n"
     "\\bibitem{Sjo85} T.~Sj\\\"ostrand, Phys.~Lett.~{\\bf 157B} (1985) 321;\n"
     "\\bibitem{Ben85} M.~Bengtsson, T.~Sj\\\"ostrand and M.~van~Zijl, "
     "Z.~Phys.~{\\bf C32} (1986) 67\n"
     "\\bibitem{Miu99} G.~Miu and T.~Sj\\\"ostrand, "
     "Phys.~Lett.~{\\bf B449} (1999) 313");

  static Reference<ShowerHandler,SpaceShowerHandler>
    interfaceSpaceShowerHandler
    ("SpaceShower",
     "Pointer to a Pythia7::SpaceShowerHandler object used to perform the "
     "showering of a space-like parton.",
     &ShowerHandler::theSpaceShower, true, false, true, true, false);

  static Reference<ShowerHandler,TimeShowerHandler> interfaceTimeShowerHandler
    ("TimeShower",
     "Pointer to a Pythia7::TimeShowerHandler object used to perform the "
     "showering of a time-like parton.",
     &ShowerHandler::theTimeShower, true, false, true, true, false);

  static Switch<ShowerHandler,bool> interfaceFSROnISR
    ("FSROnISR",
     "If space-like showers are switched on, should the partons produced "
     "in these undergo time-like cascade?",
     &ShowerHandler::addFSROnISR, true, true, false);
  static SwitchOption interfaceFSROnISRYes
    (interfaceFSROnISR,
     "Yes",
     "Time-like showers are added",
     true);
  static SwitchOption interfaceFSROnISRNo
    (interfaceFSROnISR,
     "No",
     "Time-like showers are not added.",
     false);

  static Switch<ShowerHandler,bool> interfaceFullHistory
    ("FullHistory",
     "Should all branchings be added to the event record?",
     &ShowerHandler::theFullHistory, false, true, false);
  static SwitchOption interfaceFullHistoryNo
    (interfaceFullHistory,
     "No",
     "Only final partons are added without information about the full "
     "branching.",
     false);
  static SwitchOption interfaceFullHistoryYes
    (interfaceFullHistory,
     "Yes",
     "All branchings are added.",
     true);

  static Parameter<ShowerHandler,long> interfaceMaxTries
    ("MaxTries",
     "The maximum number of attempts to cascade a sub-process before giving "
     "up and throwing an eventerror.",
     &ShowerHandler::maxTries, 1000, 1, 0,
     true, false, Interface::lowerlim);

  interfaceTimeShowerHandler.rank(10);
  interfaceSpaceShowerHandler.rank(9);

}

