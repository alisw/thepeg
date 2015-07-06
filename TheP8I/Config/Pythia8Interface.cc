// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Pythia8Interface class.
//
#include "Pythia8Interface.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/StepHandler.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

using namespace TheP8I;

Pythia8Interface::Pythia8Interface(): pythia(NULL), hooks(NULL), doHooks(false) {}

Pythia8Interface::~Pythia8Interface() {
  if ( pythia ) delete pythia;
  if ( hooks ) delete hooks;
}

bool Pythia8Interface::initialized = false;

string Pythia8Interface::installDir = PYTHIA8_DIR;

RndmEngine Pythia8Interface::rnd;

void Pythia8Interface::Init() {
  if ( initialized ) return;

  if ( installDir[installDir.length() - 1] != '/' ) installDir += '/';

  initialized = true;

}

void Pythia8Interface::errorlist() {
  pythia->info.errorStatistics();
}

bool Pythia8Interface::go() {
  try {
    return pythia->next();
  } catch ( ... ) {
      Throw<Pythia8ExecException>()
	<< "Pythia8 threw an exception while executing 'Pythia::next()'."
	<< Exception::runerror;
  }
  return false;
}

void Pythia8Interface::
setParameters(const Interfaced & handler, const vector<string> & additional) {
  if ( !pythia ) return;
  InterfaceMap ifs = BaseRepository::getInterfaces(typeid(handler));
  for ( InterfaceMap::iterator it = ifs.begin(); it != ifs.end(); ++it ) {
    string name = it->first;
    string::size_type i = name.find('_');
    ostringstream cmd;
    if ( i == string::npos ) continue;
    while ( ( i = name.find('_') ) != string::npos ) name[i] = ':';
    if ( const SwitchBase * si = dynamic_cast<const SwitchBase *>(it->second) ) {
      if ( si->get(handler) == si->def(handler) ) continue;
      cmd << name << " = " << si->get(handler);
    }
    else if ( const ParameterTBase<double> * pi =
	      dynamic_cast<const ParameterTBase<double> *>(it->second) ) {
      if ( pi->tget(handler) == pi->tdef(handler) ) continue;
      cmd << name << " = " << pi->tget(handler);
    }
    else if ( const ParameterTBase<int> * pi =
	      dynamic_cast<const ParameterTBase<int> *>(it->second) ) {
      if ( pi->tget(handler) == pi->tdef(handler) ) continue;
      cmd << name << " = " << pi->tget(handler);
    }
    else
      continue;
    pythia->readString(cmd.str());
  }
  for ( int i = 0, N = additional.size(); i < N; ++i )
    pythia->readString(additional[i]);
}

void Pythia8Interface::debug() {
  event().list();
}

void Pythia8Interface::
init(const Interfaced & handler, const vector<string> & additional) {
  if ( pythia ) delete pythia;
  if( hooks ) delete hooks;

  Init();

  CurrentGenerator::Redirect stdout(cout);

  pythia = new Pythia8::Pythia(installDir + "xmldoc");

  pythia->setRndmEnginePtr(&rnd);

  CurrentGenerator eg;
  for ( ParticleMap::const_iterator it = eg.current().particles().begin();
	it != eg.current().particles().end(); ++it ) {
    ParticleData & pd = *it->second;
    ostringstream oss;
    if ( pd.id() == 0 ) 	Throw<Pythia8InitException>()
	  << "The particle " << pd.name()
	  << " was found to have an PDG id of zero." << Exception::warning;
    oss << pd.id();
    if ( pythia->particleData.isParticle(pd.id()) )
      oss << " all ";
    else
      oss << " new ";
    oss << pd.PDGName() << " ";
    if ( pd.CC() ) {
      if ( pd.id() < 0 ) continue;
      ParticleData & apd = *pd.CC();
      oss << apd.PDGName() << " ";
      if ( pd.iSpin() != apd.iSpin() ||
	   pd.iColour() != -apd.iColour() ||
	   pd.iCharge() != -apd.iCharge() ||
	   pd.mass() != apd.mass() ||
	   pd.width() != apd.width() ||
	   pd.massMax() != apd.massMax() ||
	   pd.massMin() != apd.massMin() ||
	   pd.cTau() != apd.cTau() )
	Throw<Pythia8InitException>()
	  << "The particle " << pd.PDGName()
	  << " did not have the the same properties as its anti-particle. "
	  << "Pythia8 will be initialized with the anti-particle having "
	  << "the same properties as the particle." << Exception::warning;
    }
    else {
      oss << "void ";
    }
    oss << pd.iSpin() << " "
	<< pd.iCharge() << " "
	<< ( pd.iColour() == PDT::Colour3? 1:
	     ( pd.iColour() == PDT::Colour3bar? -1:
	       ( pd.iColour() == PDT::Colour8? 2: 0) ) ) << " "
	<< pd.mass()/GeV << " "
	<< pd.width()/GeV << " "
	<< pd.massMin()/GeV << " "
	<< pd.massMax()/GeV << " "
	<< pd.cTau()/millimeter;
    if ( !pythia->particleData.readString(oss.str(), false) )
      Throw<Pythia8InitException>()
	<< "Something went wrong when initializing the particle data table "
	<< "of Pythia8. The settings line was:\n" << oss.str()
	<< "\n" << Exception::runerror;
    pythia->particleData.hasChanged(pd.id(), false);
  }

  setParameters(handler, additional);
  
  if(doHooks){
    hooks = new RopeUserHooks();
    pythia->setUserHooksPtr(hooks);
  }

  pythia->init();

}


int Pythia8Interface::addParticle(tPPtr p, int status, int mother1, int mother2) {
  if ( particleIndex.included(p) ) return particleIndex(p);
  int idx = event().size();
  particleIndex(idx, p);
  int pos = event().append(p->id(), status, mother1, mother2, 0, 0,		
			   addColourLine(p->colourLine()),
			   addColourLine(p->antiColourLine()),
			   p->momentum().x()/GeV, p->momentum().y()/GeV,
			   p->momentum().z()/GeV, p->momentum().e()/GeV,
			   p->momentum().mass()/GeV);
  if ( !(p->vertex() == LorentzPoint()) ) {
    event()[pos].vProd(p->vertex().x()/millimeter, p->vertex().y()/millimeter,
		       p->vertex().z()/millimeter, p->vertex().t()/millimeter);
    event()[pos].tau(p->lifeLength().tau()/millimeter);
  }

  return idx;
}

int Pythia8Interface::addColourLine(tColinePtr c) {
  if ( colourIndex.included(c) ) return colourIndex(c);
  int ctag = event().nextColTag();
  colourIndex(ctag, c);
  return ctag;
}

PPtr Pythia8Interface::getParticle(int idx) {
  if ( particleIndex.included(idx) ) return particleIndex.find(idx);
  tcPDPtr pd = CurrentGenerator::current().getParticleData(event()[idx].id());
  if ( !pd ) return PPtr();
  PPtr p = pd->produceParticle(Lorentz5Momentum(event()[idx].px()*GeV,
						event()[idx].py()*GeV,
						event()[idx].pz()*GeV,
						event()[idx].e()*GeV,
						event()[idx].m()*GeV));
  if ( event()[idx].hasVertex() )
    p->setVertex(LorentzPoint(event()[idx].xProd()*millimeter,
			      event()[idx].yProd()*millimeter,
			      event()[idx].zProd()*millimeter,
			      event()[idx].tProd()*millimeter));
  if ( event()[idx].tau() > 0.0 && event()[idx].m() > 0.0 )
    p->setLifeLength(event()[idx].tau()*millimeter
		     *p->momentum()/p->momentum().mass());
    
			      
  tColinePtr c = colourIndex(event()[idx].col());
  if ( c ) c->addColoured(p);
  c = colourIndex(event()[idx].acol());
  if ( c ) c->addAntiColoured(p);
  particleIndex(idx, p);
  return p;

}
