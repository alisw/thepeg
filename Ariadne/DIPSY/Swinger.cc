// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Swinger class.
//

#include "Swinger.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/Parameter.h"
#include "CPUTimer.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Swinger.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

Swinger::~Swinger() {}

void Swinger::testGenerate(set<DipolePtr> & dips,double miny, double maxy) {
  cout << "Testing Swinger....." << endl;
  set<DipolePtr>::iterator it;
  DipoleState testState = DipoleState();
  for(it = dips.begin();it!=dips.end();it++) {
    testState.addDipole(*(*it));
    (*it)->dipoleState(& testState);
  }
  testState.sortDipoles();
  thestMode = true;
  ofstream youtput("SwingerTest.dat");
  int ysteps = 100;
  double ystepsize = 0.1;
  vector<int> ybin(ysteps,0);
  double y;
  it=dips.begin();
  for(int i=0;i<1000000;i++) {
    (*it)->generatedY(maxy);
    generate(*(*it),miny, maxy,true);
    if((*it)->swingDipole()) {
      y = (*it)->generatedY();
      if(y-miny<ysteps*ystepsize)
        ybin[int((y-miny)/ystepsize)]++;
    }
    (*it)->reset();
  }
  for(int i=0;i<ysteps;i++) {
    youtput << (i+0.5)*ystepsize << '\t' << ybin[i] << endl; 
  }
  youtput.close();
  thestMode = false;
  cout << "Done testing Swinger." 
       << endl;

}

bool Swinger::forceGenerate(Dipole & dipole, double ymax) const {
  const vector<tDipolePtr> & candidates =
    dipole.dipoleState().swingCandidates(dipole.colour());
  for (vector<tDipolePtr>::size_type i = 0; i < candidates.size(); i++ ) {
    if ( candidates[i] == & dipole ) continue;
    if ( candidates[i]->neighbors().first == & dipole || 
       candidates[i]->neighbors().second == & dipole ) continue;		
    InvEnergy2 a = dipole.partons().first->dist2(*dipole.partons().second);
    InvEnergy2 b = candidates[i]->partons().first->dist2(*(candidates[i]->partons().second));
    InvEnergy2 c = dipole.partons().second->dist2(*candidates[i]->partons().first);
    InvEnergy2 d = dipole.partons().first->dist2(*candidates[i]->partons().second);
    a = sqr(Current<DipoleEventHandler>()->rMax())*
        sqr(exp(sqrt(a)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    b = sqr(Current<DipoleEventHandler>()->rMax())*
        sqr(exp(sqrt(b)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    c = sqr(Current<DipoleEventHandler>()->rMax())*
        sqr(exp(sqrt(c)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    d = sqr(Current<DipoleEventHandler>()->rMax())*
        sqr(exp(sqrt(d)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    double yi =  -log( UseRandom::rnd() )*c*d/(a*b);
    if ( yi < dipole.generatedY() || !(dipole.swingDipole()) ) {
      if( yi < ymax || ymax == 0.0 ) { //Don't do ridiculously improbable swings.
        dipole.swingDipole(candidates[i]);
        dipole.generatedY(yi);
      }
     }
  }
  if( (dipole.swingDipole()) ) return true;
  else return false;
}

void Swinger::generateFS(Dipole & dipole, double miny, double maxy) const {
	//old implementation
  bool force = !dipole.hasGen();
  const vector<tDipolePtr> & candidates =
    dipole.dipoleState().swingCandidates(dipole.colour());
  vector<tDipolePtr>::size_type i = 0;
  while ( i < candidates.size() && candidates[i++] != & dipole );
  if ( i == candidates.size() && candidates[i - 1] != & dipole ) 
    Throw<SwingConsistencyException>()
      << "FSswinger not found among candidates" << Exception::abortnow;
  for (; i < candidates.size(); ++i ) {
    if(candidates[i]->children().first)
      cout << "OH NOES, parents among the swing candidates!!!! O_o" << endl;
    if ( dipole.neighbors().first == candidates[i] ||
	 dipole.neighbors().second == candidates[i] )
      continue;
    if ( !force && candidates[i]->hasGen() ) continue;
    // Generate a y
    double R = -log( UseRandom::rnd() );
    double amp = swingAmpFS(dipole.partons(),
     			    candidates[i]->partons(), miny/GeV2);
    //    double amp = swingAmpFS(dipole, *candidates[i], miny/GeV2);
    double yi = Constants::MaxRapidity;
    if ( miny*amp + R < Constants::MaxRapidity*amp ) yi = miny + R/amp;
    if ( yi < dipole.generatedY() || !dipole.hasGen() ) {
      dipole.swingDipole(candidates[i]);
      dipole.generatedY(yi);
      dipole.recoilSwing(false);
     }
  }
}

InvEnergy2 Swinger::
swingDistanceFS(const Parton & p1, const Parton & p2, InvEnergy2 time) const {
  //  static CPUClock cpuclock("DIPSY::Swinger::swingDinstanceFS");
  //  CPUTimer timer(cpuclock);
  static DebugItem notime("DIPSY::NoSwingTime", 60);
  if ( notime ) time = ZERO;
  Energy2 E2 = p1.plus()*p2.minus() + p1.minus()*p2.plus();
  Parton::Point x1 = p1.position() - p2.position() + time*p1.pT();
  Parton::Point x2 = p1.position() - p2.position() - time*p2.pT();
  InvEnergy2 xxScalar = x1.x()*x2.x() + x1.y()*x2.y();
  Energy2 ppScalar = p1.pT().x()*p2.pT().x() + p1.pT().y()*p2.pT().y();
  InvEnergy2 a = sqr(time)*E2 + xxScalar - sqr(time)*ppScalar;
  if ( a <= ZERO ) return ZERO;
  double asa = Current<DipoleEventHandler>()->alphaS(sqrt(a));
  InvEnergy rmax = Current<DipoleEventHandler>()->rMax();
  return sqr(rmax)/asa*sqr(exp(sqrt(a)/rmax) - 1.0);
}

double Swinger::
swingAmpFS(InvEnergy2 a, InvEnergy2 b, InvEnergy2 c, InvEnergy2 d) const {
  static double eps = 0.0000000000001;
  if ( a <= ZERO  || b <= ZERO || c <= ZERO || d <= ZERO ) return eps;
  if ( sqrt(a*b) < sqrt(c)*sqrt(d)*sqrt(theLambda*eps) ) return eps;
  return a*b/(c*d)*theLambda;
}


double Swinger::swingAmpFS(const pair<tPartonPtr, tPartonPtr> firstDip,
			   const pair<tPartonPtr, tPartonPtr> secondDip,
			   InvEnergy2 time) const {
  //parton start in their original transverse position and z = 0, then move at
  //lightspeed in a direction determined by the momentum.
  static double eps = 0.0000000000001;
  tPartonPtr p11 = firstDip.first;
  tPartonPtr p12 = firstDip.second;
  tPartonPtr p21 = secondDip.first;
  tPartonPtr p22 = secondDip.second;

  //these squared distances a,b,c,d between partons goes as (x1+t*pt1) - (x2+t*pt2) for small t,
  //and as invariant mass times t for large t
  Energy2 E2 = p11->plus()*p12->minus() + p11->minus()*p12->plus();
  Parton::Point x1 = p11->position() - p12->position() + time*p11->pT();
  Parton::Point x2 = p11->position() - p12->position() - time*p12->pT();
  InvEnergy2 xxScalar = x1.x()*x2.x() + x1.y()*x2.y();
  Energy2 ppScalar = p11->pT().x()*p12->pT().x() + p11->pT().y()*p12->pT().y();
  InvEnergy2 a = sqr(time)*E2 + xxScalar - sqr(time)*ppScalar;

  E2 = p21->plus()*p22->minus() + p21->minus()*p22->plus();
  x1 = p21->position() - p22->position() + time*p21->pT();
  x2 = p21->position() - p22->position() - time*p22->pT();
  xxScalar = x1.x()*x2.x() + x1.y()*x2.y();
  ppScalar = p21->pT().x()*p22->pT().x() + p21->pT().y()*p22->pT().y();
  InvEnergy2 b = sqr(time)*E2 + xxScalar - sqr(time)*ppScalar;

  E2 = p21->plus()*p12->minus() + p21->minus()*p12->plus();
  x1 = p21->position() - p12->position() + time*p21->pT();
  x2 = p21->position() - p12->position() - time*p12->pT();
  xxScalar = x1.x()*x2.x() + x1.y()*x2.y();
  ppScalar = p21->pT().x()*p12->pT().x() + p21->pT().y()*p12->pT().y();
  InvEnergy2 c = sqr(time)*E2 + xxScalar - sqr(time)*ppScalar;

  E2 = p11->plus()*p22->minus() + p11->minus()*p22->plus();
  x1 = p11->position() - p22->position() + time*p11->pT();
  x2 = p11->position() - p22->position() - time*p22->pT();
  xxScalar = x1.x()*x2.x() + x1.y()*x2.y();
  ppScalar = p11->pT().x()*p22->pT().x() + p11->pT().y()*p22->pT().y();
  InvEnergy2 d = sqr(time)*E2 + xxScalar - sqr(time)*ppScalar;

  if ( a <= ZERO  || b <= ZERO || c <= ZERO || d <= ZERO ) return eps;

  //normal confinement correction
  InvEnergy rmax = Current<DipoleEventHandler>()->rMax();
  double asa = Current<DipoleEventHandler>()->alphaS(sqrt(a));
  double asb = Current<DipoleEventHandler>()->alphaS(sqrt(b));
  double asc = Current<DipoleEventHandler>()->alphaS(sqrt(c));
  double asd = Current<DipoleEventHandler>()->alphaS(sqrt(d));
  a = sqr(rmax)/asa*sqr(exp(sqrt(a)/rmax) - 1.0);
  b = sqr(rmax)/asb*sqr(exp(sqrt(b)/rmax) - 1.0);
  c = sqr(rmax)/asc*sqr(exp(sqrt(c)/rmax) - 1.0);
  d = sqr(rmax)/asd*sqr(exp(sqrt(d)/rmax) - 1.0);

  //normal amplitude
  if ( sqrt(a*b) < sqrt(c)*sqrt(d)*sqrt(theLambda*eps) ) return eps;
  return a*b/(c*d)*theLambda;

}

void Swinger::generate(Dipole & dipole, double miny, double maxy, bool force) const {
  double yi = max(max(dipole.partons().first->y(),
		      dipole.partons().second->y()), miny);
  if ( dipole.partons().second == dipole.partons().first )
    Throw<ColourIndexException>()
      << "Found inconsistent colour indices in DIPSY. Event discarded."
      << Exception::eventerror;
  if( dipole.children().first)
    Throw<SwingConsistencyException>()
      << "OH NOES, parents among the swing candidates!!!! O_o"
      << Exception::abortnow;
  const vector<tDipolePtr> & candidates =
    dipole.dipoleState().swingCandidates(dipole.colour(), !force);
  vector<tDipolePtr>::const_iterator it = candidates.begin();
  if ( force ) while ( it != candidates.end() && *it++ != &dipole );
  for (; it != candidates.end(); ++it ) {
    if ( !force && (**it).hasGen() ) continue;
    double yn = yi - log(UseRandom::rnd())/swingAmp(dipole, **it);
    if ( yn < dipole.generatedY() ) {
      dipole.swingDipole(*it);
      dipole.generatedY(yn);
      dipole.recoilSwing(false);
    }
  }
}

double Swinger::swingAmp(const Dipole & firstDip,
			 const Dipole & secondDip) const {
  InvEnergy rmax = Current<DipoleEventHandler>()->rMax();
  InvEnergy2 a = firstDip.swingCache;
  if ( a < ZERO ) {
    a = firstDip.partons().first->dist2(*firstDip.partons().second);
    a = sqr(rmax)/(Current<DipoleEventHandler>()->alphaS(sqrt(a)))*
      sqr(exp(sqrt(a)/rmax) - 1.0);
    firstDip.swingCache = a;
  }
  InvEnergy2 b = secondDip.swingCache;
  if ( b < ZERO ) {
    b = secondDip.partons().first->dist2(*(secondDip.partons().second));
    b = sqr(rmax)/(Current<DipoleEventHandler>()->alphaS(sqrt(b)))*
      sqr(exp(sqrt(b)/rmax) - 1.0);
    secondDip.swingCache = b;
  }
  InvEnergy2 c = firstDip.partons().second->dist2(*secondDip.partons().first);
  c = sqr(rmax)/(Current<DipoleEventHandler>()->alphaS(sqrt(c)))*
    sqr(exp(sqrt(c)/rmax) - 1.0);
  InvEnergy2 d = firstDip.partons().first->dist2(*secondDip.partons().second);
  d = sqr(rmax)/(Current<DipoleEventHandler>()->alphaS(sqrt(d)))*
    sqr(exp(sqrt(d)/rmax) - 1.0);
  return a*b/(c*d)*theLambda;
}

double Swinger::swingAmp(const pair<tPartonPtr, tPartonPtr> firstDip,
			  const pair<tPartonPtr, tPartonPtr> secondDip) const {
  InvEnergy rmax = Current<DipoleEventHandler>()->rMax();
  InvEnergy2 a = firstDip.first->dist2(*firstDip.second);
  InvEnergy2 b = secondDip.first->dist2(*(secondDip.second));
  InvEnergy2 c = firstDip.second->dist2(*secondDip.first);
  InvEnergy2 d = firstDip.first->dist2(*secondDip.second);
  a = sqr(rmax)/
    (Current<DipoleEventHandler>()->alphaS(sqrt(a)))*
    sqr(exp(sqrt(a)/rmax) - 1.0);
  b = sqr(rmax)/
    (Current<DipoleEventHandler>()->alphaS(sqrt(b)))*
    sqr(exp(sqrt(b)/rmax) - 1.0);
  c = sqr(rmax)/
    (Current<DipoleEventHandler>()->alphaS(sqrt(c)))*
    sqr(exp(sqrt(c)/rmax) - 1.0);
  d = sqr(rmax)/
    (Current<DipoleEventHandler>()->alphaS(sqrt(d)))*
    sqr(exp(sqrt(d)/rmax) - 1.0);
  return a*b/(c*d)*theLambda;

}

bool Swinger::checkMaxY(DipolePtr d1, DipolePtr d2, double maxy) const {
  tPartonPtr p1 = d1->partons().first;
  tPartonPtr p2 = d1->partons().second;
  tPartonPtr p3 = d2->partons().first;
  tPartonPtr p4 = d2->partons().second;
  TransverseMomentum rec12 = (p1->position() - p2->position())/
    (p1->position() - p2->position()).pt2();
  TransverseMomentum rec34 = (p3->position() - p4->position())/
    (p3->position() - p4->position()).pt2();
  TransverseMomentum rec14 = (p1->position() - p4->position())/
    (p1->position() - p4->position()).pt2();
  TransverseMomentum rec32 = (p3->position() - p2->position())/
    (p3->position() - p2->position()).pt2();
  if ( log( (p1->pT() - rec12 + rec14).pt()/p1->plus()) > maxy )
    return false;
  if ( log( (p2->pT() + rec12 - rec32).pt()/p2->plus()) > maxy )
    return false;
  if ( log( (p3->pT() - rec34 + rec32).pt()/p3->plus()) > maxy )
    return false;
  if ( log( (p4->pT() + rec34 - rec14).pt()/p4->plus()) > maxy )
    return false;
  return true;
}

void Swinger::recombineFS(Dipole & d1) const {
  Dipole & d2 = *d1.swingDipole();

  pair<tPartonPtr, tPartonPtr> d1p = d1.partons();
  pair<tPartonPtr, tPartonPtr> d2p = d2.partons();
  pair<tDipolePtr, tDipolePtr> d1n = d1.neighbors();
  pair<tDipolePtr, tDipolePtr> d2n = d2.neighbors();

  d1.partons(make_pair(d1p.first, d2p.second));
  d1.neighbors(make_pair(d1n.first, d2n.second));
  if ( d1n.first ) d1n.first->secondNeighbor(&d1);
  if ( d2n.second ) d2n.second->firstNeighbor(&d1);

  d2.partons(make_pair(d2p.first, d1p.second));
  d2.neighbors(make_pair(d2n.first, d1n.second));
  if ( d2n.first ) d2n.first->secondNeighbor(&d2);
  if ( d1n.second ) d1n.second->firstNeighbor(&d2);

  d1.partons().first->dipoles(make_pair(d1.neighbors().first, &d1));
  d1.partons().second->dipoles(make_pair(&d1, d1.neighbors().second));
  d2.partons().first->dipoles(make_pair(d2.neighbors().first, &d2));
  d2.partons().second->dipoles(make_pair(&d2, d2.neighbors().second));

  d1.reset();
  d2.reset();
  d1.touch();
  d2.touch();

}

void Swinger::recombine(Dipole & d1) const {
  DipoleState & state = d1.dipoleState();
  Dipole & d2 = *d1.swingDipole();
  
  if(d2.children().second) {
    Throw<SwingConsistencyException>()
      << "swinging with someone thats emitted already (second child). "
      << "(partons at " << d1.partons().first->y() <<  ", "
      << d1.partons().second->y() << "; " << d2.partons().first->y()
      <<  ", " << d2.partons().second->y() << ")" << Exception::eventerror;
    state.diagnosis(true);
  }
  else if(d2.children().first) {
    Throw<SwingConsistencyException>()
      << "swinging with someone thats emitted already (first child). "
      << "(partons at " << d1.partons().first->y() <<  ", "
      << d1.partons().second->y() << "; " << d2.partons().first->y()
      <<  ", " << d2.partons().second->y() << ")" << Exception::eventerror;
    state.diagnosis(true);
  }

  DipolePtr d11 = d1.children().first = state.createDipole();
  DipolePtr d22 = d2.children().first = state.createDipole();
  
  d11->partons(make_pair(d1.partons().first, d2.partons().second));
  d11->neighbors(make_pair(d1.neighbors().first, d2.neighbors().second));
  if ( d1.neighbors().first ) d1.neighbors().first->secondNeighbor(d11);
  if ( d2.neighbors().second ) d2.neighbors().second->firstNeighbor(d11);
  
  d22->partons(make_pair(d2.partons().first, d1.partons().second));
  d22->neighbors(make_pair(d2.neighbors().first, d1.neighbors().second));
  if ( d2.neighbors().first ) d2.neighbors().first->secondNeighbor(d22);
  if ( d1.neighbors().second ) d1.neighbors().second->firstNeighbor(d22);

  d11->colour(d1.colour());
  d22->colour(d2.colour());

  d11->partons().first->dipoles(make_pair(d11->neighbors().first,d11));
  d11->partons().second->dipoles(make_pair(d11,d11->neighbors().second));
  d22->partons().first->dipoles(make_pair(d22->neighbors().first,d22));
  d22->partons().second->dipoles(make_pair(d22,d22->neighbors().second)); 

  if (d1.DGLAPsafe() || d2.DGLAPsafe() ) {
    d11->DGLAPsafe(true);
    d22->DGLAPsafe(true);
  }

  if ( d1.interacted() ) d11->interacted(d1.interacted());
  if ( d2.interacted() ) d22->interacted(d2.interacted());

  d1.reset();
  d2.reset();          //are these really needed?
}

void Swinger::persistentOutput(PersistentOStream & os) const {
  os << ounit(theLambda, 1.0);
}

void Swinger::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theLambda, 1.0);
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<Swinger,HandlerBase>
  describeDIPSYSwinger("DIPSY::Swinger", "libAriadne5.so libDIPSY.so");

void Swinger::Init() {

  static ClassDocumentation<Swinger> documentation
    ("There is no documentation for the Swinger class");

  static Parameter<Swinger,double> interfaceLambda
    ("Lambda",
     "The frequency of the swings.",
     &Swinger::theLambda, 1.0, 1.0, 0.0, 0.0,
     true, false, Interface::lowerlim);

}

