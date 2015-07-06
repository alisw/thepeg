// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleSwinger class.
//

#include "FSGluonEmission.h"
#include "DipoleSwinger.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/StateDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/EnumIO.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "../../DIPSY/CPUTimer.h"

using namespace Ariadne5;

DipoleSwinger::DipoleSwinger():
  lambda(1.0), Rmax(3.5*InvGeV*hbarc), linear(false), sizeOpt(0),
  theRhoCut(-1.0*GeV), theMaxRho(-2.0*GeV), stringMassCut(ZERO), disallowedrho(0.0*GeV) {}

DipoleSwinger::~DipoleSwinger() {}

IBPtr DipoleSwinger::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleSwinger::fullclone() const {
  return new_ptr(*this);
}

void DipoleSwinger::TauRatio::plot(Time dt1, Time dt2, int N) const {
  double rmax = maximum(dt1, dt2);
  for ( int i = 0; i < N; ++i ) {
    Time dt = dt1 + (i + 0.5)*(dt2 - dt1)/double(N);
    cerr << dt/femtometer << '\t'
	 << rat(dt)/rmax << '\t'
	 << ounit(sqrt(d1tau2(dt)), femtometer) << '\t'
	 << ounit(sqrt(d2tau2(dt)), femtometer) << '\t'
	 << ounit(sqrt(r1tau2(dt)), femtometer) << '\t'
	 << ounit(sqrt(r2tau2(dt)), femtometer) << endl;
  }
}

void DipoleSwinger::TauRatio::debug() const {
  cerr << sqrt(-d1tau2.dx2)/femtometer << '\t' <<  sqrt(-d2tau2.dx2)/femtometer << endl;
  cerr << sqrt(-r1tau2.dx2)/femtometer << '\t' <<  sqrt(-r2tau2.dx2)/femtometer << endl;
  cerr << "   inf " << rat(ZERO) << '\t' << chris(ZERO) << endl;
  cerr << "100.00 " << rat(hbarc/(100.0*GeV)) << '\t' << chris(hbarc/(100.0*GeV)) << endl;
  cerr << " 10.00 " << rat(hbarc/(10.0*GeV)) << '\t' << chris(hbarc/(10.0*GeV)) << endl;
  cerr << "  1.00 " << rat(hbarc/(1.0*GeV)) << '\t' << chris(hbarc/(1.0*GeV)) << endl;
  cerr << "  0.10 " << rat(hbarc/(0.1*GeV)) << '\t' << chris(hbarc/(0.1*GeV)) << endl;
  cerr << "  0.01 " << rat(hbarc/(0.01*GeV)) << '\t' << chris(hbarc/(0.01*GeV)) << endl;
  cerr << "  0.00 " << rat(hbarc/(0.000001*GeV)) << '\t' << chris(hbarc/(0.000001*GeV)) << endl;
  cerr << "       " << srat << endl;
  cerr << sqrt(s12)/GeV << '\t' <<  sqrt(s34)/GeV << endl;
  cerr << sqrt(s14)/GeV << '\t' <<  sqrt(s32)/GeV << endl;
}

DipoleSwinger::TauRatio::TauRatio(const QCDDipole & d1, const QCDDipole & d2,
				  Length Rmaxin, int optin)
  : d1tau2(d1.iPart(), d1.oPart(), optin, Rmaxin),
    d2tau2(d2.iPart(), d2.oPart(), optin, Rmaxin),
    r1tau2(d1.iPart(), d2.oPart(), optin, Rmaxin),
    r2tau2(d2.iPart(), d1.oPart(), optin, Rmaxin),
    Rmax(Rmaxin), sizeOpt(optin) {
  s12 = (d1.iPart()->momentum() + d1.oPart()->momentum()).m2();
  s34 = (d2.iPart()->momentum() + d2.oPart()->momentum()).m2();
  s14 = (d1.iPart()->momentum() + d2.oPart()->momentum()).m2();
  s32 = (d2.iPart()->momentum() + d1.oPart()->momentum()).m2();
  srat = (s12*s34)/(s14*s32);
  // LorentzMomentum ptot = d1.iPart()->momentum() + d1.oPart()->momentum() +
  //                        d2.iPart()->momentum() + d2.oPart()->momentum();
}

double DipoleSwinger::TauRatio::rat(Time dt) const {
  if ( Rmax == ZERO )
    return d1tau2(dt)*d2tau2(dt)/(r1tau2(dt)*r2tau2(dt));
  if ( sizeOpt == 0 )
    return sqr((exp(sqrt(d1tau2(dt))/abs(Rmax)) - 1.0)*(exp(sqrt(d2tau2(dt))/abs(Rmax)) - 1.0)/
	       ((exp(sqrt(r1tau2(dt))/abs(Rmax)) - 1.0)*(exp(sqrt(r2tau2(dt))/abs(Rmax)) - 1.0)));
  else if ( sizeOpt == 2 )
    return exp(-sqr(min(0.0*meter, d1tau2.db() + d2tau2.db() - r1tau2.db() - r2tau2.db())/dt))*
      (d1tau2(dt)*d2tau2(dt))/(r1tau2(dt)*r2tau2(dt));    
  else if ( Rmax > ZERO )
    return exp((d1tau2.db2() + d2tau2.db2() - r1tau2.db2() - r2tau2.db2())/sqr(Rmax))*
      (d1tau2(dt)*d2tau2(dt))/(r1tau2(dt)*r2tau2(dt));
  else
    return exp(-(d1tau2.db() + d2tau2.db() - r1tau2.db() - r2tau2.db())/Rmax)*
      (d1tau2(dt)*d2tau2(dt))/(r1tau2(dt)*r2tau2(dt));
}

double DipoleSwinger::TauRatio::maximum(Time dt1, Time dt2) const {
  // return d1tau2.minmax(dt1, dt2).second*d2tau2.minmax(dt1, dt2).second/
  // 	(r1tau2.minmax(dt1, dt2).first*r2tau2.minmax(dt1, dt2).first);
  double rmax = max(rat(dt1), rat(dt2));
  if ( sizeOpt >= 1 ) return 1.01*rmax;
  double rmax2 = max(rat(r1tau2.extreme(dt1, dt2)), rat(r2tau2.extreme(dt1, dt2)));
  if ( rmax2 > rmax ) rmax = rmax2;
  return 1.1*rmax;
}

bool DipoleSwinger::canHandle(const DipoleBase & e) const {
  tcStateDipPtr d = dynamic_ptr_cast<tcStateDipPtr>(&e);
  if ( !d ) return false;
  return true;
}

bool DipoleSwinger::overrides(const EmitterBase &, DipoleBase &) const {
  return false;
}

bool DipoleSwinger::touched(const DipoleBase & dipole) const {
  return dipole.state()->touched() || dipole.touched() ||
    dipole.state() != lastState || dipole.state()->uniqueId != lastStateId;
}

Energy DipoleSwinger::rhoCut() const {
  return theRhoCut > ZERO? theRhoCut: EmitterBase::rhoCut();
}

EmPtr DipoleSwinger::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  static DebugItem logme("Ariadne5::DipoleSwinger", 6);
  static DebugItem swinglast("Ariadne5::SwingLast", 60);
  static CPUClock cpuclock("Ariadne5::DipoleSwinger::generate");
  static CPUClock cpuclock1("Ariadne5::DipoleSwinger::generate1");
  static CPUClock cpuclock2("Ariadne5::DipoleSwinger::generate2");
  static CPUClock cpuclock3("Ariadne5::DipoleSwinger::generate3");
  CPUTimer timer(cpuclock);


  // Only for debugging.
  if ( swinglast && rhomax > 1.0*GeV ) return EmPtr();

  // Some swings results in too small strings and are rejected, these
  // should be disallowed in the following generations. All such
  // dipole pairs are collected in disallowed set until an accepted
  // emission (one which lowers the rhomax) is found.
  if ( rhomax != disallowedrho ) disallowed.clear();
  if ( lastSelected && lastSelected->state == Emission::reverted ) {
    disallowed.insert(lastSelected->dipoles);
    disallowedrho = rhomax;
  }

  if ( maxRho() > rhoCut() ) rhomax = min(rhomax, maxRho());
  //  if ( rhomax <= rhomin ) return EmPtr();

  // Here we collect all dipoes which may swing in different sets
  // depending on their colour index. Touched dipoles are alse stored
  // in a separate set.
  const set<tDBPtr> & activeset = dipole.state()->activeDipoles();
  map<ColourIndex, vector<tQCDPtr> > cactive;
  map<ColourIndex, vector<tQCDPtr> > ctouched;

  for ( set<tDBPtr>::const_iterator it = activeset.begin();
	it != activeset.end(); ++it )
    if ( tQCDPtr dp = dynamic_ptr_cast<tQCDPtr>(*it) ) {
      cactive[dp->colourIndex()].push_back(dp);
      if ( dp->touched() || dp->iPart()->touched() || dp->oPart()->touched() )
	ctouched[dp->colourIndex()].push_back(dp);
    }

  CPUTimer timer1(cpuclock1);

  // Collect initial values in time (inverse scale) variables.
  Time dtmax = hbarc/rhomin;
  Time dtmin = hbarc/rhomax;
  Time dtmaxcut = hbarc/rhoCut();
  Time t0 = hbarc/Current<AriadneHandler>()->pTCut();

  // Now check how much we need to redo since last time.
  bool redoFull = false;

  // If the previous emission was a forced swing of a tiny string, we
  // need to redo everything, because everything was not considered
  // last time.
  if ( lastSelected && lastSelected->state == Emission::performed && lastSelected->forced )
    redoFull = true;

  // If this is a new dipole state we clearly have to redo everything.
  if ( redoFull || dipole.state() != lastState || dipole.state()->uniqueId != lastStateId ) {
    redoFull = true;
    lastState = dipole.state();
    lastStateId = dipole.state()->uniqueId;
    cache.clear();
    findTinyStrings();
    lastSelected = DipSwPtr();
  }

  // If the state dipole was touched it means the lastSelected swing was performed.
  if ( redoFull || dipole.touched() ) lastSelected = DipSwPtr();

  // If either dipole in the previous generation was touched can't use the last generated swing.
  if ( lastSelected && ( lastSelected->dipoles.first->touched() ||
			 lastSelected->dipoles.second->touched() ) )
    lastSelected = DipSwPtr();

  // Just for debugging.
  if ( logme && redoFull )
    generator()->log() << "DipoleSwinger: Starting new event." << endl;
  if ( logme ) generator()->log()
		 << "DipoleSwinger: Starting new swing search." << endl;

  bool forceTiny = !tinyString.empty();
  for ( map<ColourIndex, vector<tQCDPtr> >::iterator ciit = cactive.begin();
	ciit != cactive.end(); ++ciit ) {
    const vector<tQCDPtr> & active = ciit->second;
    if ( active.empty() ) continue;
    const vector<tQCDPtr> & touched = ctouched[active[0]->colourIndex()];
    int ti1 = touched.empty()? -1: 0;
    
    for ( int i1 = 0, N = active.size(); i1 < N; ++i1  ) {
      CPUTimer timer2(cpuclock2);
      if ( ti1 >= 0 && touched[ti1] == active[i1] ) ++ti1;
      QCDDipole & d1 = *active[i1];
      bool inTiny = false;
      if ( forceTiny ) inTiny = member(tinyString, &d1);

      tcQCDPtr nextnext = tcQCDPtr();
      if ( stringMassCut > ZERO && d1.next() && d1.next()->sdip() < sqr(stringMassCut) )
	nextnext = d1.next()->next();
      tcQCDPtr prevprev = tcQCDPtr();
      if ( stringMassCut > ZERO && d1.prev() && d1.prev()->sdip() < sqr(stringMassCut) )
	prevprev = d1.prev()->prev();

      bool redo = redoFull || d1.touched() ||
	d1.iPart()->touched() || d1.oPart()->touched();

      pair<CacheMap::iterator,bool> cit = cache.insert(make_pair(&d1, DipSwPtr()));
      if ( redo ) cit.first->second = DipSwPtr();
      DipSwPtr sel = cit.first->second;
      if ( !sel || ( sel->dipoles.second->touched() ||
		     sel->dipoles.second->iPart()->touched() ||
		     sel->dipoles.second->oPart()->touched() ||
		     sel->dipoles.second->colourIndex() != d1.colourIndex() ||
		     sel->dipoles.second == nextnext ||
		     sel->dipoles.second == prevprev ) ) {
	redo = true;
	sel = cit.first->second = DipSwPtr();
      }
      if ( !redo && sel && disallowed.find(sel->dipoles) != disallowed.end() ) {
	redo = true;
	sel = cit.first->second = DipSwPtr();
      }

      const vector<tQCDPtr> * secondp = &active;
      int i2 = i1 + 1;
      if ( !redo ) {
	if ( ti1 < 0 ) {
	  select(sel);
	  continue;
	}
	secondp = &touched;
	i2 = ti1;
      }
      const vector<tQCDPtr> & second = *secondp;

      CPUTimer timer3(cpuclock3);
      for ( int N2 = second.size(); i2 < N2; ++i2 ) {
	if ( second[i2] == nextnext || second[i2] == prevprev ) continue;
	if ( disallowed.find(make_pair(active[i1], second[i2])) != disallowed.end() ) continue;
	QCDDipole & d2 = *second[i2];
	if ( d1.iPart() == d2.oPart() || d2.iPart() == d1.oPart() ) {
	  cerr << "Bullocks! Try to swing a gluon into colour singlet!" << endl;
	}

	// If we are forcing a swing on a tiny string we should only
	// consider swings where one dipole is in the tyny string and
	// the other is not.
	if ( forceTiny && inTiny == member(tinyString, &d2) ) continue;

	TauRatio tauRatio(d1, d2, Rmax, sizeOpt);
	// LorentzMomentum ptot = d1.iPart()->momentum() + d1.oPart()->momentum() +
	//                        d2.iPart()->momentum() + d2.oPart()->momentum();
	double ratmax = tauRatio.maximum(dtmin, dtmax);
	double Cinv = -1.0/(lambda*ratmax);
	Time dt = dtmin;
	Time dtcut = dtmax;
	if ( sel ) dtcut = min(dtcut, hbarc/sel->rho);
 
	do {
	  double logR = log(UseRandom::rnd());
	  if ( linear )
	    dt += logR*Cinv*t0;
	  else {
	    if ( logR*Cinv > log(2.0*dtmaxcut/dt) ) dt = dtmaxcut*2.0;
	    else dt *= exp(logR*Cinv);
	  }
	  if ( dt >= dtcut ) break;
	  if ( tauRatio(dt) > ratmax ) {
	    cerr << "Ooops failed to overestimate swing probability: "
		 << tauRatio(dt)/ratmax << endl;
	  }
	} while ( dt < dtcut && tauRatio(dt) < UseRandom::rnd()*ratmax );


	if ( !sel ) sel = new_ptr(DipoleSwing(*this, dipole, d1, d2, dt));
	else if ( dt < dtcut ) sel->setup(d1, d2, dt);

      }

      select(cit.first->second = sel);

    }
  }

  

  if ( lastSelected && forceTiny ) {
    lastSelected->rho = rhomax;
    lastSelected->forced = true;
  }
  if ( lastSelected && lastSelected->rho < rhomin ) lastSelected = DipSwPtr();
    
  return lastSelected;

}

void DipoleSwinger::findTinyStrings() const {
  tinyString.clear();
  if ( stringMassCut <= ZERO ) return;
  MinCmp<Energy> cut(stringMassCut);
  set<tDBPtr> tdipoles = lastState->activeDipoles();
  set<tcDBPtr> dipoles(tdipoles.begin(), tdipoles.end());

  while ( !dipoles.empty() ) {
    set<tcDBPtr> currentstring;
    if ( tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(*dipoles.begin()) ) {
      dipoles.erase(d);
      pair<tcQCDPtr,tcQCDPtr> range = DipoleState::StringEnds(d);
      tParPtr p1 = range.first->iPart();
      tParPtr p2 = range.second->oPart();
      LorentzMomentum psum;
      while ( range.first ) {
	currentstring.insert(range.first);
	psum += range.first->iPart()->momentum();
	if ( range.first == range.second )
	  range.first = tcQCDPtr();
	else
	  range.first = range.first->next();
      }
      if ( p1 != p2 ) psum += p2->momentum();
      if ( cut(psum.m()) ) {
	tinyString = currentstring;
      }
    } else {
      dipoles.erase(dipoles.begin());
    }
  }
}      
  


void DipoleSwinger::swing(tQCDPtr d1, tQCDPtr d2) {
  // Swing the o-partons.
  tParPtr tmp = d1->oPart();
  d1->oPart(d2->oPart());
  d2->oPart(tmp);
  // The new o-partons inhreits the previous colour line.
  tColinePtr cltmp = d1->oPart()->origICol();
  d1->oPart()->origICol(d2->oPart()->origICol());
  d2->oPart()->origICol(cltmp);
  // Fix pointers to neighboring dipoles.
  tQCDPtr dtmp = d1->next();
  d1->next(d2->next());
  d2->next(dtmp);
  if ( d1->next() ) d1->next()->prev(d1);
  if ( d2->next() ) d2->next()->prev(d2);
}

bool DipoleSwinger::
perform(const Emission & emission) const {
  static DebugItem stat("Ariadne5::DipoleSwinger::Statistics", 6);
  static double sum = 0.0;
  static CPUClock cpuclock("Ariadne5::DipoleSwinger::perform");
  CPUTimer timer(cpuclock);
  
  const DipoleSwing & e = dynamic_cast<const DipoleSwing &>(emission);

  if ( e.dipoles.first->colourIndex() != e.dipoles.second->colourIndex() )
    cerr << "Ooops, tried to swing dipoles of diferent colours!" << endl;

  double lam0 = 0.0;
  Energy2 s12 = ZERO;
  Energy2 s34 = ZERO;

  if ( stat ) {
    lam0 = log(e.dipoles.first->sdip()/GeV2) + log(e.dipoles.second->sdip()/GeV2);
    TauRatio tauRatio(*e.dipoles.first, *e.dipoles.second, Rmax, sizeOpt);
    cerr << "ratio " << tauRatio(hbarc/e.rho);
    s12 = e.dipoles.first->sdip();
    s34 = e.dipoles.second->sdip();
  }


  swing(e.dipoles.first, e.dipoles.second);

  if ( stringMassCut > ZERO && !e.forced ) {
    if ( Utilities::sumMomentum(e.dipoles.first->string()).m() < stringMassCut ) return false;
    if ( Utilities::sumMomentum(e.dipoles.second->string()).m() < stringMassCut ) return false;
  }

  if ( stat ) {
    double lam = log(e.dipoles.first->sdip()/GeV2) + log(e.dipoles.second->sdip()/GeV2);
    cerr << " (" << s12*s34/(e.dipoles.first->sdip()*e.dipoles.second->sdip())
	 << "), swing diff " << lam - lam0
	 << " (total " << (sum += lam - lam0) << ")" << endl;
  }

  e.dipoles.first->touch();
  e.dipoles.second->touch();
  e.dipole->touch();

  return true;

}


void DipoleSwinger::revert(const Emission & emission) const {
  static CPUClock cpuclock("Ariadne5::DipoleSwinger::revert");
  CPUTimer timer(cpuclock);
  const DipoleSwing & e = dynamic_cast<const DipoleSwing &>(emission);

  swing(e.dipoles.first, e.dipoles.second);
  e.dipoles.second->untouch();
  e.dipole->untouch();

}




// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleSwinger::persistentOutput(PersistentOStream & os) const {
  os << lambda << ounit(Rmax, femtometer) << oenum(linear) << sizeOpt
     << ounit(theRhoCut, GeV) << ounit(theMaxRho, GeV) << ounit(stringMassCut, GeV)
     << tinyString << lastSelected;
}

void DipoleSwinger::persistentInput(PersistentIStream & is, int) {
  is >> lambda >> iunit(Rmax, femtometer) >> ienum(linear) >> sizeOpt
     >> iunit(theRhoCut, GeV) >> iunit(theMaxRho, GeV) >> iunit(stringMassCut, GeV)
     >> tinyString >> lastSelected;
  cache.clear();
  lastState = tDipoleStatePtr();
  lastStateId = 0;
}

string DipoleSwinger::setRmax(string cmd) {
  istringstream is(cmd);
  double x;
  is >> x;
  Rmax = x*InvGeV*hbarc;
  return "";
}

DescribeClass<DipoleSwinger,EmitterBase>
describeAriadne5DipoleSwinger("Ariadne5::DipoleSwinger", "libAriadne5.so");

void DipoleSwinger::Init() {

  static ClassDocumentation<DipoleSwinger> documentation
    ("The DipoleSwinger class implements the final-state swing method"
     "for colour reconnections.");

  static Parameter<DipoleSwinger,double> interfaceLambda
    ("Lambda",
     "The frequency of the swings.",
     &DipoleSwinger::lambda, 1.0, 1.0, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Parameter<DipoleSwinger,Length> interfaceRmax
    ("Rmax",
     "The typical hadronic size used to regularize large dipoles in the swing.",
     &DipoleSwinger::Rmax, femtometer, 3.5*InvGeV*hbarc, 0.0*femtometer, 0.0*femtometer,
     true, false, Interface::nolimits);

  static Command<DipoleSwinger> interfaceSetRmax
    ("SetRmax",
     "Set the value of <interface>Rmax</interface> with a value given in units of inverse GeV.",
     &DipoleSwinger::setRmax, true);


  static Switch<DipoleSwinger,bool> interfaceLinear
    ("Linear",
     "Normally we use logarithmic evolution in time. With this switch we can instead use linear evolution.",
     &DipoleSwinger::linear, false, true, false);
  static SwitchOption interfaceLinearLogarithmic
    (interfaceLinear,
     "Logarithmic",
     "Use logarithmic evolution in time.",
     false);
  static SwitchOption interfaceLinearLinear
    (interfaceLinear,
     "Linear",
     "Use linear evolution in time.",
     true);

  static Switch<DipoleSwinger,int> interfaceSizeOpt
    ("SizeOpt",
     "Different options for how to define dipole size is the swing.",
     &DipoleSwinger::sizeOpt, 0, true, false);
  static SwitchOption interfaceSizeOptNaive
    (interfaceSizeOpt,
     "Naive",
     "Assuming a global time given by the inverse transverse momentum, use the "
     "Lorentz-invariant actual distance.",
     0);
  static SwitchOption interfaceSizeOptInvariantMass
    (interfaceSizeOpt,
     "InvariantMass",
     "Use only the invariant mass of the dipoles.",
     1);
  static SwitchOption interfaceSizeOptBSupInvariantMass
    (interfaceSizeOpt,
     "BSupInvariantMass",
     "Use only the invariant mass of the dipoles but suppress swinging to large dipoles.",
     2);
  static SwitchOption interfaceSizeOptInvariantMassR
    (interfaceSizeOpt,
     "NaiveR",
     "As for Naive but Rmax limit only applied to original vertex positions.",
     -1);

  static Parameter<DipoleSwinger,Energy> interfaceRhoCut
    ("RhoCut",
     "Cutoff in the evolution variable. If zero or less, use the global QCD shower cutoff.",
     &DipoleSwinger::theRhoCut, GeV, -1.0*GeV, -1.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<DipoleSwinger,Energy> interfaceRhoMax
    ("MaxRho",
     "The maximum allowed value of the evolution variable, Used to disallow perturbative "
     "swings. Inactive if below <interface>RhoCut</interface>.",
     &DipoleSwinger::theMaxRho, GeV, -2.0*GeV, -2.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<DipoleSwinger,Energy> interfaceStringMassCut
    ("StringMassCut",
     "Cutoff against too small string. No swing is allowed to produce strings with invarint "
     "mass below this value.",
     &DipoleSwinger::stringMassCut, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

}

