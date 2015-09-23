// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleXSec class.
//

#include "DipoleXSec.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "Parton.h"
#include "DipoleEventHandler.h"
#include "RealParton.h"
#include "RealPartonState.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DipoleXSec.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "gsl/gsl_sf_bessel.h"

using namespace ThePEG;
using namespace DIPSY;

DipoleXSec::~DipoleXSec() {}

DipoleInteraction::DipoleInteraction(Dipole & dlin, Dipole & drin,
				     const ImpactParameters &  bin, int ordering)
  : dips(&dlin, &drin), dnext(dlin.neighbors().first, drin.neighbors().first),
    dprev(dlin.neighbors().second, drin.neighbors().second),
    b(&bin), sints(tSPartonPtr(), tSPartonPtr()),
    f2(0.0), uf2(0.0), norec(false), kt(ZERO), id(0), status(UNKNOWN),
    intOrdering(ordering) {

 // If different colours interaction is either between c-c or
 // cbar-cbar. If the same colour we only haver c-cbar.
  ints = make_pair(dlin.partons().first, drin.partons().first);
  spec = make_pair(dlin.partons().second, drin.partons().second);
  if ( dlin.colour() == drin.colour() ) {
    //    if ( UseRandom::rnd() > 0.5  )
    if ( UseRandom::rndbool()  )
      swap(ints.second, spec.second);
    else
      swap(ints.first, spec.first);
  }      
  if ( Current<DipoleEventHandler>()->eventFiller().compat ) {
    // *** TODO *** debugging to be removed!
    InvEnergy2 rr11 = bin.dist2(*dlin.partons().first, *drin.partons().first);
    InvEnergy2 rr22 = bin.dist2(*dlin.partons().second, *drin.partons().second);
    double dummy = (rr11 + rr22)*GeV2;
    while ( dummy > 100 ) dummy /= 10;
    if ( 1000000.0*dummy - long(1000000.0*dummy) > 0.5 ) {
      swap(ints.first, spec.first);
      swap(ints.second, spec.second);
    }
  } else {
    //    if ( UseRandom::rnd() < 0.5 ) {
    if ( UseRandom::rndbool() ) {
      swap(ints.first, spec.first);
      swap(ints.second, spec.second);
    }
  }


  Parton::Point r = bin.difference(ints.first->position(), ints.second->position()); 
  d2 = min(min(r.pt2(), dlin.size2()), drin.size2());
  //  kt = 1.0/sqrt(d2);
  rec =  -Current<DipoleEventHandler>()->emitter().pTScale()*r/r.pt2();
  kt = rec.pt();
  if ( ints.first->shadow() ) {
    sints = make_pair(ints.first->shadow()->resolveInteraction(d2, spec.first),
		      ints.second->shadow()->resolveInteraction(d2, spec.second));
  }
}

void DipoleInteraction::prepare() const {
  sints = make_pair(ints.first->shadow()->resolveInteraction(d2, spec.first),
		    ints.second->shadow()->resolveInteraction(d2, spec.second));
  sints.first->insertInteraction(id);
  sints.second->insertInteraction(id);
}

void DipoleInteraction::debug() const {
  sints.first->memememe = sints.second->memememe = true;
  TransverseMomentum rrec;
  if ( !norec ) rrec = rec;
  cerr << setw(15) << rrec.x()/GeV << setw(15) << rrec.y()/GeV << endl;
  dips.first->dipoleState().debugShadowTree();
  rrec = b->invRotatePT(-rrec);
  cerr << setw(15) << rrec.x()/GeV << setw(15) << rrec.y()/GeV << endl;
  dips.second->dipoleState().debugShadowTree();
  sints.first->memememe = sints.second->memememe = false;
  list<PartonPtr> plist = dips.first->dipoleState().getPartons();
  LorentzMomentum sum;
  for ( list<PartonPtr>::iterator it = plist.begin(); it != plist.end(); ++it )
    if ( (**it).onShell() ) sum += (**it).momentum();
  cerr << "left:  "
       << setw(15) << sum.x()/GeV << setw(15) << sum.y()/GeV
       << setw(15) << sum.z()/GeV << setw(15) << sum.t()/GeV << endl;
  plist = dips.second->dipoleState().getPartons();
  sum = LorentzMomentum();
  for ( list<PartonPtr>::iterator it = plist.begin(); it != plist.end(); ++it )
    if ( (**it).onShell() ) sum += (**it).momentum();
  TransverseMomentum sumt = b->rotatePT(TransverseMomentum(sum.x(), sum.y()));
  cerr << "right:  "
       << setw(15) << sumt.x()/GeV << setw(15) << sumt.y()/GeV
       << setw(15) << sum.z()/GeV << setw(15) << sum.t()/GeV << endl;
}

void DipoleInteraction::fail(int i ) const {
  int off = i < 0? 4: 0;
  if ( i > 0 ) status = ORDERING;
  else {
    if ( ints.first == dips.first->partons().first &&
	 ints.second == dips.second->partons().first )
      i = 10;
    else if ( ints.first == dips.first->partons().first )
      i = 11;
    else if ( ints.second == dips.second->partons().first )
      i = 12;
    else
      i = 13;
  }
  i += off;
  
  ++ofail[i];
  if ( uf2 > 0.2 ) ++o1fail[i];
}

// *** TODO *** Remove status, it is only for debugging.
DipoleInteraction::Status DipoleInteraction::check(int mode) const {
  static DebugItem breakme("DIPSY::Stop", 6);
  if ( mode >= 0 && breakme ) breakThePEG();
  status = UNKNOWN;
  sints = make_pair(ints.first->shadow()->resolveInteraction(d2, spec.first),
		    ints.second->shadow()->resolveInteraction(d2, spec.second));
  ShadowParton::Propagator ppl =
    sints.first->propagator(d2, spec.first, mode);
  ShadowParton::Propagator ppr =
    sints.second->propagator(d2, spec.second, mode);

  fail(-1);
  if ( ppl.fail || ppr.fail ) return status = PROPFAIL;
  LorentzMomentum pl = ppl.p;
  LorentzMomentum pr = ppr.p;
  if ( mode >= 0 ) {
    TransverseMomentum ptr0(pr.x(), pr.y());
    checkShadowMomentum(pl + lightCone(pr.minus(), pr.plus(),
				       b->rotatePT(ptr0)));
  }
  TransverseMomentum rrec;
  if ( !norec ) rrec = rec;
  TransverseMomentum ptl =
    TransverseMomentum(pl.x(), pl.y()) + rrec;
  Energy2 mtl2 = ptl.pt2() + sqr(sints.first->mass());
  TransverseMomentum ptr =
    TransverseMomentum(pr.x(), pr.y()) + b->invRotatePT(-rrec);
  Energy2 mtr2 = ptr.pt2() + sqr(sints.second->mass());
  Energy PP = pl.plus() + pr.minus();
  Energy PM = pl.minus() + pr.plus();
  if ( PP <= ZERO || PP*PM + mtl2 - mtr2 < ZERO ||
       PM <= ZERO || PP*PM + mtr2 - mtl2 < ZERO ) return status = KINEFAIL;
  Energy2 sqrl = (sqr(PP*PM + mtl2 - mtr2) - 4.0*mtl2*PM*PP)/sqr(PM);
  Energy2 sqrr = (sqr(PP*PM + mtr2 - mtl2) - 4.0*mtr2*PM*PP)/sqr(PP);
  if ( sqrl < ZERO || sqrr < ZERO ) return status = KINEFAIL;
  sints.first->pTplus(ptl, 0.5*(PP + mtl2/PM - mtr2/PM + sqrt(sqrl)));
  sints.second->pTplus(ptr, 0.5*(PM + mtr2/PP - mtl2/PP + sqrt(sqrr)));

  if ( intOrdering <= 2 ) {
    if ( PP > min(ppl.colpos, ppl.acopos) ||
	 PM > min(ppr.colpos, ppr.acopos) ) return status = ORDERING;
  } else if ( intOrdering == 4 ) {
    if ( orderfail(ppl, ppr) ) return status = ORDERING;
  }
  if ( !id ) return status = ACCEPTED;

  sints.first->setOnShell(mode);
  sints.second->setOnShell(mode);
  if ( mode >= 0 ) {
    sints.first->interacting(id);
    sints.first->acceptInteraction(id);
    sints.second->interacting(id);
    sints.second->acceptInteraction(id);
  }
  if ( mode > 0 ) {
    sints.first->original()->interact(true);
    sints.second->original()->interact(true);
  }

  if ( mode >= 0 ) checkShadowMomentum();

  return status = ACCEPTED;
}

void DipoleInteraction::checkShadowMomentum(const LorentzMomentum & pin) const {
  static DebugItem checkkinematics("DIPSY::CheckKinematics", 6);
  if ( !checkkinematics ) return;
  Sum20Momentum sum20;
  dips.first->dipoleState().checkShadowMomentum(sum20);
  dips.second->dipoleState().checkShadowMomentum(sum20, b);
  sum20 -= pin;

  if ( !sum20 ) {
    Throw<InteractionKinematicException>()
      << "Shadow trees had inconsistent momentum" << Exception::warning;
    debug();
  }
}

bool DipoleInteraction::orderfail(const ShadowParton::Propagator & ppl,
				  const ShadowParton::Propagator & ppr) const {
  pair<bool,bool> parcol(ints.first == dips.first->partons().first,
			 ints.second == dips.second->partons().first);
  Energy2 ptl = ppl.p.perp2();
  Energy2 ptr = ppr.p.perp2();
  Energy2 ptm = (ppl.p - sints.first->momentum()).perp2();
  Energy2 ptlm = max(ptl, ptm);
  Energy2 ptrm = max(ptr, ptm);

  fail(0);
 
  if ( parcol.first &&  disorder(cending(ppl.acoptp, ptl, ptm),
				 ppl.acopos, ppl.aconeg, sints.first) )
  // If left is coloured it should be ordered wrt. previous emissions
  // on the anti-colour line.
    fail(1);
  if ( !parcol.first && disorder(cending(ppl.colptp, ptl, ptm),
				 ppl.colpos, ppl.colneg, sints.first) )
  // And vice vers.
    fail(2);

  if ( parcol.second &&  disorder(cending(ppr.acoptp, ptr, ptm),
				  ppr.acopos, ppr.aconeg, sints.second) )
  // If right is coloured it should be ordered wrt. previous emissions
  // on the anti-colour line.
    fail(3);
  if ( !parcol.second && disorder(cending(ppr.colptp, ptr, ptm),
				  ppr.colpos, ppr.colneg, sints.second) )
  // And vice vers.
    fail(4);

  if ( parcol.first != parcol.second &&  disorder(cending(ptl, ptm, ptr)) )
  // If colour--anto-colour combinations,left and right should be
  // ordered wrt. eachother
    fail(5);

  if ( parcol.first && parcol.second ) {
    // If both partons are coloured they should be ordered wrt. the
    // others incoming emissions on the colour line.
    if ( disorder(cending(ppl.colptp, ptlm, ppr.colptp),
		  sints.first, ppr.colneg, ppr.colpos) )
      fail(6);
    if ( disorder(cending(ppr.colptp, ptrm, ppl.colptp),
		  sints.second, ppl.colneg, ppl.colpos) )
      fail(7);
  }
  if ( !parcol.first && !parcol.second ) {
    // ... and vice versa for anti.colour - anti-colour
    if ( disorder(cending(ppl.acoptp, ptlm, ppr.acoptp),
		  sints.first, ppr.aconeg, ppr.acopos) )
      fail(8);
    if ( disorder(cending(ppr.acoptp, ptrm, ppl.acoptp),
		  sints.second, ppl.aconeg, ppl.acopos) )
      fail(9);
  }

  return ( status == ORDERING );
}

void DipoleInteraction::reject() const {
  ints.first->shadow()->rejectInteraction(id);
  ints.second->shadow()->rejectInteraction(id);
}

void DipoleInteraction::accept() const {
  sints = make_pair(ints.first->shadow()->resolveInteraction(d2, spec.first),
		    ints.second->shadow()->resolveInteraction(d2, spec.second));
}

InvEnergy DipoleXSec::rMax() const {
  return theRMax > 0.0*InvGeV? theRMax: Current<DipoleEventHandler>()->rMax();
}

double DipoleXSec::fSinFn(const Parton::Point & rho1, const Parton::Point & rho2,
			  const TransverseMomentum & pt) const {
  
  if ( sinFunction == 1 ) {
    double r1 = rho1.pt2()*pt.pt2();
    double r2 = rho2.pt2()*pt.pt2();
    return r1*r2/(4.0*(r1 + 1.0)*(r2 + 1.0));
  }

  double s1 = pt.x()*rho1.x() + pt.y()*rho1.y();
  double s2 = pt.x()*rho2.x() + pt.y()*rho2.y();
  return sqr(sin(s1)*sin(s2));

}

double DipoleXSec::
fij(const pair<tPartonPtr, tPartonPtr> left,
    const pair<tPartonPtr, tPartonPtr> right,
    const ImpactParameters & b, bool veto) const {
  Nfij++;

  tcPartonPtr p11 = left.first;
  tcPartonPtr p12 = left.second;
  tcPartonPtr p21 = right.first;
  tcPartonPtr p22 = right.second;

  //TODO: keep only interaction 0, as that is the only one supported anyway
  if ( theInteraction == 0 || theInteraction == 1 || theInteraction == 3 || theInteraction == 4 ) {
    InvEnergy2 rr11 = b.dist2(*p11, *p21);
    InvEnergy2 rr21 = b.dist2(*p12, *p21);
    InvEnergy2 rr12 = b.dist2(*p11, *p22);
    InvEnergy2 rr22 = b.dist2(*p12, *p22);

    if ( veto && kinematicsVeto(left, right, b) ) {
      return 0.0;
    }

    TransverseMomentum pt;
    double pTScale = Current<DipoleEventHandler>()->emitter().pTScale();

    InvEnergy rscale = sqrt(min(min(min(rr12, rr21),min(rr11, rr22)),
				min(p11->dist2(*p12), p21->dist2(*p22))));

    double fudgeME= 1.0;

    if ( theInteraction == 0 ) {  
      pair<bool, bool> ints;
      ints = int0Partons(p11, p12, p21, p22, b);

      Parton::Point r1 = b.difference((ints.first ? p11->position():p12->position()),
				      (ints.second ? p21->position():p22->position()));

      if ( r1.pt2() < min(p11->dist2(*p12), p21->dist2(*p22)) &&
	   Current<DipoleEventHandler>()->fudgeME() ) {
	double deltay = ( ints.first? p11->y(): p12->y() ) +
	                ( ints.second? p21->y(): p22->y() );
	fudgeME = 1.0 - 1.0/(1.0 + cosh(deltay));
	fudgeME *= Current<DipoleEventHandler>()->fudgeFactorME();
      }

      pt = pTScale*r1/r1.pt2();
      rscale = sqrt(min(min(r1.pt2(),p11->dist2(*p12)), p21->dist2(*p22)));
    }

    if ( theInteraction == 1 ) {  //4 parton distance
      Parton::Point r1 = b.difference(p11->position(), p21->position());
      Parton::Point r2 = b.difference(p11->position(), p22->position());
      Parton::Point r3 = b.difference(p12->position(), p21->position());
      Parton::Point r4 = b.difference(p12->position(), p22->position());
      pt = pTScale*(r1/sqr(r1.pt()) + r2/sqr(r2.pt()) +
				       r3/sqr(r3.pt()) + r4/sqr(r4.pt()))*0.35;
    }

    if ( theInteraction == 3 ) {  //2 parton distance
      Parton::Point r1 = b.difference(p11->position(), p22->position());
      Parton::Point r2 = b.difference(p12->position(), p21->position());
      pt = pTScale*(r1/sqr(r1.pt()) + r2/sqr(r2.pt()))*0.6;
    }

    Parton::Point rho1 = (p12->position() - p11->position())/2.0;
    Parton::Point rho2 = b.rotate(p22->position() - p21->position())/2.0;

    return 8.0*fudgeME*sqr(Current<DipoleEventHandler>()->alphaS(rscale))*
      fSinFn(rho1, rho2, pt)*
      sqr(sqr(pt.pt()))/sqr(sqr(pt.pt()) + sqr(pTScale/rMax()));
  }

  return 0.0;
}

double DipoleXSec::
fij(const Dipole & dleft, const Dipole & dright,
    const ImpactParameters & b, bool veto) const {
  pair<tPartonPtr, tPartonPtr> left = dleft.partons();
  pair<tPartonPtr, tPartonPtr> right = dright.partons();
  pair<bool, bool> ints;

  Nfij++;

  double fudgeME = 1.0;
  tcPartonPtr p11 = left.first;
  tcPartonPtr p12 = left.second;
  tcPartonPtr p21 = right.first;
  tcPartonPtr p22 = right.second;
  ints = int0Partons(p11, p12, p21, p22, b);
  //TODO: keep only interaction 0, as that is the only one supported anyway
  if ( theInteraction == 0 || theInteraction == 1 || theInteraction == 3 || theInteraction == 4 ) {
    InvEnergy2 rr11 = b.dist2(*p11, *p21);
    InvEnergy2 rr21 = b.dist2(*p12, *p21);
    InvEnergy2 rr12 = b.dist2(*p11, *p22);
    InvEnergy2 rr22 = b.dist2(*p12, *p22);

    if ( veto && kinematicsVeto(dleft, dright, b, ints) ) {
      return 0.0;
    }

    TransverseMomentum pt;
    double pTScale = Current<DipoleEventHandler>()->emitter().pTScale();

    InvEnergy rscale = sqrt(min(min(min(rr12, rr21),min(rr11, rr22)),
				min(p11->dist2(*p12), p21->dist2(*p22))));

    if ( theInteraction == 0 ) {  

      Parton::Point r1 = b.difference((ints.first ? p11->position():p12->position()),
				      (ints.second ? p21->position():p22->position()));

      if ( r1.pt2() < min(p11->dist2(*p12), p21->dist2(*p22)) &&
	   Current<DipoleEventHandler>()->fudgeME() ) {
	double deltay = ( ints.first? p11->y(): p12->y() ) +
	                ( ints.second? p21->y(): p22->y() );
	fudgeME = 1.0 - 1.0/(1.0 + cosh(deltay));
      }

      if ( Current<DipoleEventHandler>()->fudgeME() > 1 ) {
	if ( dleft.colour() == dright.colour() ) fudgeME *= 9.0/2.0;
	else fudgeME *= 9.0/16.0;
      }

      pt = pTScale*r1/sqr(r1.pt());
      rscale = sqrt(min(min(sqr(r1.pt()),p11->dist2(*p12)), p21->dist2(*p22)));
    }

    if ( theInteraction == 1 ) {  //4 parton distance
      Parton::Point r1 = b.difference(p11->position(), p21->position());
      Parton::Point r2 = b.difference(p11->position(), p22->position());
      Parton::Point r3 = b.difference(p12->position(), p21->position());
      Parton::Point r4 = b.difference(p12->position(), p22->position());
      pt = pTScale*(r1/sqr(r1.pt()) + r2/sqr(r2.pt()) +
				       r3/sqr(r3.pt()) + r4/sqr(r4.pt()))*0.35;
    }

    if ( theInteraction == 3 ) {  //2 parton distance
      Parton::Point r1 = b.difference(p11->position(), p22->position());
      Parton::Point r2 = b.difference(p12->position(), p21->position());
      pt = pTScale*(r1/sqr(r1.pt()) + r2/sqr(r2.pt()))*0.6;
    }

    Parton::Point rho1 = (p12->position() - p11->position())/2.0;
    Parton::Point rho2 = b.rotate(p22->position() - p21->position())/2.0;

    return 8.0*fudgeME*sqr(Current<DipoleEventHandler>()->alphaS(rscale))*
      fSinFn(rho1, rho2, pt)*
      sqr(sqr(pt.pt()))/sqr(sqr(pt.pt()) + sqr(pTScale/rMax()));
  }

  return 0.0;
}

DipoleInteraction DipoleXSec::
fij(const ImpactParameters & b, Dipole & dleft, Dipole & dright,
    bool veto) const {
  DipoleInteraction di(dleft, dright, b, theIntOrdering);

  Nfij++;

  double fudgeME = 1.0;

  if ( veto && kinematicsVeto(di) ) return di;

  Energy m0 = Current<DipoleEventHandler>()->emitter().pTScale()/rMax();


  Parton::Point r1 = b.difference(di.ints.first->position(), di.ints.second->position());

  if ( r1.pt2() < min(dleft.size2(), dright.size2()) &&
       Current<DipoleEventHandler>()->fudgeME() ) {
    double deltay =  di.ints.first->y() + di.ints.second->y();
    fudgeME = 1.0 - 1.0/(1.0 + cosh(deltay));
    fudgeME *= Current<DipoleEventHandler>()->fudgeFactorME();
  }
  if ( Current<DipoleEventHandler>()->fudgeME() > 1 ) {
    if ( dleft.colour() == dright.colour() ) fudgeME *= 9.0/2.0;
    else fudgeME *= 9.0/16.0;
  }

  Parton::Point rho1 = di.dips.first->vSize()/2.0;
  Parton::Point rho2 = b.rotate(di.dips.second->vSize())/2.0;

  di.f2 = 16.0*fudgeME*sqr(Current<DipoleEventHandler>()->alphaS(sqrt(di.d2)))*
    fSinFn(rho1, rho2, di.rec)*sqr(di.rec.pt2())/sqr(di.rec.pt2() + sqr(m0));
  di.uf2 = unitarize(di.f2);

  return di;
}

bool DipoleXSec::kinematicsVeto(const pair<tPartonPtr, tPartonPtr> left,
				const pair<tPartonPtr, tPartonPtr> right,
				const ImpactParameters &  b) const {

  pair<pair<bool, bool>, pair<bool, bool> > ints = doesInt(left, right, b);
  if ( Current<DipoleEventHandler>()->eventFiller().compat ) {
    pair<bool,bool> int0 =
      int0Partons(left.first, left.second, right.first, right.second, b);
    ints = make_pair(make_pair(int0.first, !int0.first),
 		     make_pair(int0.second, !int0.second));
  }

  InteractionRecoil recs = recoil(left, right, b, ints);

  if ( theInteraction == 0 ) {

    //the key partons for f_ij
    pair<bool, bool> ints0 = int0Partons(left.first, left.second, right.first, right.second, b);

    //first set up effective partons with range etc
    tPartonPtr p1 = (ints0.first ? left.first:left.second);
    tPartonPtr p2 = (ints0.second ? right.first:right.second);

    //there MAY be secondary partons in the interaction, decided by ints.
    tPartonPtr p1sec, p2sec;
    if ( ints.first.first && ints.first.second )
      p1sec = (ints0.first ? left.second:left.first);
    if ( ints.second.first && ints.second.second )
      p2sec = (ints0.second ? right.second:right.first);

    //only distance between key partons will affect range (reasonable?)
    InvEnergy range1 = sqrt(min(left.first->dist2(*left.second)/4.0, b.dist2(*p1, *p2)));
    InvEnergy range2 = sqrt(min(right.first->dist2(*right.second)/4.0, b.dist2(*p1, *p2)));
    if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 ) {
      range1 = sqrt(left.first->dist2(*left.second)/4.0);
      range2 = sqrt(right.first->dist2(*right.second)/4.0);
    }
    EffectivePartonPtr ep1 = EffectiveParton::create(*p1, range1);
    EffectivePartonPtr ep2 = EffectiveParton::create(*p2, range2);
    EffectivePartonPtr ep1sec, ep2sec;
    if ( p1sec ) ep1sec = EffectiveParton::create(*p1sec, range1);
    if ( p2sec ) ep2sec = EffectiveParton::create(*p2sec, range2);

    TransverseMomentum rec1 = (ints0.first ? recs.first.first:recs.first.second);
    TransverseMomentum rec2 = (ints0.second ? recs.second.first:recs.second.second);

    Energy pt1 = (ep1->pT() + rec1).pt();
    Energy minus1 = sqr(pt1)/ep1->plus();
    Energy pt2 = (ep2->pT() + rec2).pt();
    Energy minus2 = sqr(pt2)/ep2->plus();

    //sum up total supplied and needed LC momentum from each side
    Energy leftPlus = ep1->plus() + (p1sec ? ep1sec->plus():ZERO);
    Energy leftMinus = minus1;
    if ( p1sec ) {
    TransverseMomentum rec1sec = (ints0.first ? recs.first.second:recs.first.first);
    Energy pt1sec = (ep1sec->pT() + rec1sec).pt();
    Energy minus1sec = sqr(pt1sec)/ep1sec->plus();
    leftMinus += minus1sec;
    }
    Energy rightPlus = ep2->plus() + (p2sec ? ep2sec->plus():ZERO);
    Energy rightMinus = minus2 + (p2sec ? ep2sec->minus():ZERO);
    if ( p2sec ) {
    TransverseMomentum rec2sec = (ints0.second ? recs.second.second:recs.second.first);
    Energy pt2sec = (ep2sec->pT() + rec2sec).pt();
    Energy minus2sec = sqr(pt2sec)/ep2sec->plus();
    rightMinus += minus2sec;
    }

    if ( theIntOrdering == 2 ) {
      Energy maxRec = max(max(recs.first.first.pt() ,recs.first.second.pt()),
			  max(recs.second.first.pt() ,recs.second.second.pt()));
      if ( leftPlus*rightPlus < 16.0*sqr(maxRec) ) return true;
      else return false;
    }

    //check enough energy to set all on shell
    if ( leftPlus*rightPlus < 16.0*leftMinus*rightMinus ) {
      return true;
    }

    //take LC momentum transfers in account
    Energy plus1 = ep1->plus()*(1.0 - rightMinus/leftPlus);
    minus1 /= 1.0 - rightMinus/leftPlus;
    Energy plus2 = ep2->plus()*(1.0 - leftMinus/rightPlus);
    minus2 /= 1.0 - leftMinus/rightPlus;

    //check ordering of the key partons (secondaries can be unordered if they want)
    double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
    if ( theIntOrdering == 0 ) {
      //just check plus and minus ordering after recoils
      if ( plus1*PSInf < minus2 ) return true;
      if ( plus2*PSInf < minus1 ) return true;
    }
    else if ( theIntOrdering == 1 ) {
      //check p- ordering from both sides, as if it was evolution.
      //that is, as if no future recoils.
      double PMOrd = Current<DipoleEventHandler>()->emitter().PMinusOrdering();
      TransverseMomentum rec1 = (ints0.first ? recs.first.first:recs.first.second);
      TransverseMomentum rec2 = (ints0.second ? recs.second.first:recs.second.second);
      if ( sqr(max(rec1.pt(), rec2.pt())*PSInf) < minus1*minus2*PMOrd )
    	return true;
    }

    //check that the partons stay on their side, to avoid doublecounting
    double yInt = 0;
    if ( p1->dipoles().first ) yInt = p1->dipoles().first->dipoleState().ymax();
    else if ( p1->dipoles().second ) yInt = p1->dipoles().second->dipoleState().ymax();
    double y1 = 0.5*log(minus1/plus1);
    if ( y1 > yInt ) return true;
    double y2 = 0.5*log(plus2/minus2);
    if ( y2 < yInt ) return true;

    //should we check secondaries as well here?

    //if no veto triggered, it should not be vetoed
    return false;
  }
  else {
    cerr << "only interaction 0 is supported in kinematics veto atm, sorry" << endl;
  }
  return false;
}

bool DipoleXSec::kinematicsVeto(const Dipole & dleft,
				const Dipole & dright,
				const ImpactParameters &  b,
				const pair<bool,bool> & ints0) const {

  pair<tPartonPtr, tPartonPtr> left = dleft.partons();
  pair<tPartonPtr, tPartonPtr> right = dright.partons();
  pair<pair<bool, bool>, pair<bool, bool> > ints = doesInt(left, right, b);
  if ( Current<DipoleEventHandler>()->eventFiller().compat ) {
    ints = make_pair(make_pair(ints0.first, !ints0.first),
 		     make_pair(ints0.second, !ints0.second));
  }
  InteractionRecoil recs = recoil(left, right, b, ints, ints0);

  if ( theInteraction == 0 ) {

    //first set up effective partons with range etc
    tPartonPtr p1 = (ints0.first ? left.first:left.second);
    tPartonPtr p2 = (ints0.second ? right.first:right.second);

    //there MAY be secondary partons in the interaction, decided by ints.
    tPartonPtr p1sec, p2sec;
    if ( ints.first.first && ints.first.second )
      p1sec = (ints0.first ? left.second:left.first);
    if ( ints.second.first && ints.second.second )
      p2sec = (ints0.second ? right.second:right.first);

    //only distance between key partons will affect range (reasonable?)
    InvEnergy range1 = sqrt(min(left.first->dist2(*left.second)/4.0, b.dist2(*p1, *p2)));
    InvEnergy range2 = sqrt(min(right.first->dist2(*right.second)/4.0, b.dist2(*p1, *p2)));
    if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 ) {
      range1 = sqrt(left.first->dist2(*left.second)/4.0);
      range2 = sqrt(right.first->dist2(*right.second)/4.0);
    }
    EffectivePartonPtr ep1 = dleft.getEff(p1, range1);
    EffectivePartonPtr ep2 = dright.getEff(p2, range2);
    EffectivePartonPtr ep1sec, ep2sec;
    if ( p1sec ) ep1sec = dleft.getEff(p1sec, range1);
    if ( p2sec ) ep2sec = dright.getEff(p2sec, range2);

    TransverseMomentum rec1 = (ints0.first ? recs.first.first:recs.first.second);
    TransverseMomentum rec2 = (ints0.second ? recs.second.first:recs.second.second);

    Energy pt1 = (ep1->pT() + rec1).pt();
    Energy minus1 = sqr(pt1)/ep1->plus();
    Energy pt2 = (ep2->pT() + rec2).pt();
    Energy minus2 = sqr(pt2)/ep2->plus();

    //sum up total supplied and needed LC momentum from each side
    Energy leftPlus = ep1->plus() + (p1sec ? ep1sec->plus():ZERO);
    Energy leftMinus = minus1;
    if ( p1sec ) {
    TransverseMomentum rec1sec = (ints0.first ? recs.first.second:recs.first.first);
    Energy pt1sec = (ep1sec->pT() + rec1sec).pt();
    Energy minus1sec = sqr(pt1sec)/ep1sec->plus();
    leftMinus += minus1sec;
    }
    Energy rightPlus = ep2->plus() + (p2sec ? ep2sec->plus():ZERO);
    Energy rightMinus = minus2 + (p2sec ? ep2sec->minus():ZERO);
    if ( p2sec ) {
    TransverseMomentum rec2sec = (ints0.second ? recs.second.second:recs.second.first);
    Energy pt2sec = (ep2sec->pT() + rec2sec).pt();
    Energy minus2sec = sqr(pt2sec)/ep2sec->plus();
    rightMinus += minus2sec;
    }

    if ( theIntOrdering == 2 ) {
      Energy maxRec = max(max(recs.first.first.pt() ,recs.first.second.pt()),
			  max(recs.second.first.pt() ,recs.second.second.pt()));
      if ( leftPlus*rightPlus < 16.0*sqr(maxRec) ) return true;
      else return false;
    }

    //check enough energy to set all on shell
    if ( leftPlus*rightPlus < 16.0*leftMinus*rightMinus ) {
      return true;
    }

    //take LC momentum transfers in account
    Energy plus1 = ep1->plus()*(1.0 - rightMinus/leftPlus);
    minus1 /= 1.0 - rightMinus/leftPlus;
    Energy plus2 = ep2->plus()*(1.0 - leftMinus/rightPlus);
    minus2 /= 1.0 - leftMinus/rightPlus;

    //check ordering of the key partons (secondaries can be unordered if they want)
    double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
    if ( theIntOrdering == 0 ) {
      //just check plus and minus ordering after recoils
      if ( plus1*PSInf < minus2 ) return true;
      if ( plus2*PSInf < minus1 ) return true;
    }
    else if ( theIntOrdering == 1 ) {
      //check p- ordering from both sides, as if it was evolution.
      //that is, as if no future recoils.
      double PMOrd = Current<DipoleEventHandler>()->emitter().PMinusOrdering();
      TransverseMomentum rec1 = (ints0.first ? recs.first.first:recs.first.second);
      TransverseMomentum rec2 = (ints0.second ? recs.second.first:recs.second.second);
      if ( sqr(max(rec1.pt(), rec2.pt())*PSInf) < minus1*minus2*PMOrd )
    	return true;
    }

    //check that the partons stay on their side, to avoid doublecounting
    double yInt = 0;
    if ( p1->dipoles().first ) yInt = p1->dipoles().first->dipoleState().ymax();
    else if ( p1->dipoles().second ) yInt = p1->dipoles().second->dipoleState().ymax();
    double y1 = 0.5*log(minus1/plus1);
    if ( y1 > yInt ) return true;
    double y2 = 0.5*log(plus2/minus2);
    if ( y2 < yInt ) return true;

    //should we check secondaries as well here?

    //if no veto triggered, it should not be vetoed
    return false;
  }
  else {
    cerr << "only interaction 0 is supported in kinematics veto atm, sorry" << endl;
  }
  return false;
}

bool DipoleXSec::kinematicsVeto(const DipoleInteraction & di) const {

  if ( di.sints.first ) return di.check(-1); // Hey! We're using shadows!

  ImpactParameters b = *di.b;
  const Dipole & dleft = *di.dips.first;
  const Dipole & dright = *di.dips.second;

  //first set up effective partons with range etc
  tcPartonPtr p1 = di.ints.first;
  tcPartonPtr p2 = di.ints.second;

  //only distance between key partons will affect range (reasonable?)
  InvEnergy range1 = sqrt(min(dleft.size2()/4.0, b.dist2(*p1, *p2)));
  InvEnergy range2 = sqrt(min(dright.size2()/4.0, b.dist2(*p1, *p2)));
  EffectivePartonPtr ep1 = dleft.getEff(p1, range1);
  EffectivePartonPtr ep2 = dright.getEff(p2, range2);
  
  TransverseMomentum rec1 = di.rec;
  TransverseMomentum rec2 = di.b->invRotatePT(-di.rec);

  Energy pt1 = (ep1->pT() + rec1).pt();
  Energy minus1 = sqr(pt1)/ep1->plus();
  Energy pt2 = (ep2->pT() + rec2).pt();
  Energy minus2 = sqr(pt2)/ep2->plus();
  
    //sum up total supplied and needed LC momentum from each side
    Energy leftPlus = ep1->plus();
    Energy leftMinus = minus1;
    Energy rightPlus = ep2->plus();
    Energy rightMinus = minus2;

    if ( theIntOrdering == 2 ) {
      if ( leftPlus*rightPlus < 16.0*rec1.pt2() ) return true;
      else return false;
    }

    //check enough energy to set all on shell
    if ( leftPlus*rightPlus < 16.0*leftMinus*rightMinus ) {
      return true;
    }

    //take LC momentum transfers in account
    Energy plus1 = ep1->plus()*(1.0 - rightMinus/leftPlus);
    minus1 /= 1.0 - rightMinus/leftPlus;
    Energy plus2 = ep2->plus()*(1.0 - leftMinus/rightPlus);
    minus2 /= 1.0 - leftMinus/rightPlus;

    //check ordering of the key partons (secondaries can be unordered if they want)
    double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
    if ( theIntOrdering == 0 ) {
      //just check plus and minus ordering after recoils
      if ( plus1*PSInf < minus2 ) return true;
      if ( plus2*PSInf < minus1 ) return true;
    }

    //check that the partons stay on their side, to avoid doublecounting
    double yInt = 0;
    if ( p1->dipoles().first ) yInt = p1->dipoles().first->dipoleState().ymax();
    else if ( p1->dipoles().second ) yInt = p1->dipoles().second->dipoleState().ymax();
    double y1 = 0.5*log(minus1/plus1);
    if ( y1 > yInt ) return true;
    double y2 = 0.5*log(plus2/minus2);
    if ( y2 < yInt ) return true;

    //should we check secondaries as well here?

    //if no veto triggered, it should not be vetoed
    return false;
  
  return false;
}

double DipoleXSec::
sumf(const DipoleState & sl, const DipoleState & sr,
     const ImpactParameters & b) const {

  Nfij = 0;
  NBVeto = 0;
  NScalVeto = 0;
  scalVeto = 0.0;
  bVeto = 0.0;
  vector<tDipolePtr> dl;
  sl.extract(back_inserter(dl));
  vector<tDipolePtr> dr;
  sr.extract(back_inserter(dr));

  double sum = 0.0;
  for ( int i = 0, N = dl.size(); i < N; ++i )
    for ( int j = 0, M = dr.size(); j < M; ++j )
      sum += fij(*dl[i], *dr[j], b);
  return sum;

}

double DipoleXSec::
sumf(const ImpactParameters & b, const DipoleState & sl, const DipoleState & sr) const {

  Nfij = 0;
  NBVeto = 0;
  NScalVeto = 0;
  scalVeto = 0.0;
  bVeto = 0.0;
  vector<tDipolePtr> dl;
  sl.extract(back_inserter(dl));
  vector<tDipolePtr> dr;
  sr.extract(back_inserter(dr));

  double sum = 0.0;
  for ( int i = 0, N = dl.size(); i < N; ++i )
    for ( int j = 0, M = dr.size(); j < M; ++j ) {
      if ( (M*i + j)%8 == 0 ) UseRandom::rnd();
      DipoleInteraction di =  fij(b, *dl[i], *dr[j], false);
      if ( !di.check(-1) )
	sum += di.f2/2.0; /*** TODO: this is because we need proper f here. FIX CONFUSION */
    }
  return sum;

}

DipoleXSec::FList
DipoleXSec::flist(const DipoleState & sl, const DipoleState & sr,
		  const ImpactParameters & b) const {
  FList ret;
  vector<tDipolePtr> dl;
  sl.extract(back_inserter(dl));
  vector<tDipolePtr> dr;
  sr.extract(back_inserter(dr));
  //dont save the lowest fij. Otherwise the FList for LHC PbPb will take too much memory.
  double cutoff = min(0.00000001, 1.0/double(dl.size()*dr.size()));
  if ( dl.size()*dr.size() < 1000000 ) cutoff = 0.0000000001;
  int total = 0;
  int skipped = 0;
  for ( int i = 0, N = dl.size(); i < N; ++i ) {
    for ( int j = 0, M = dr.size(); j < M; ++j ) {
      double f = fij(*dl[i], *dr[j], b, false)*2.0;
      //extra 2 added from the non-diffractive interaction probability.
      total++;
      if ( f > cutoff &&
	   !kinematicsVeto(dl[i]->partons(), dr[j]->partons(), b) )
	ret.insert(make_pair(make_pair(f, unitarize(f)), make_pair(dl[i], dr[j])));
      else skipped++;
    }
  }
  return ret;
}

DipoleInteraction::List
DipoleXSec::flist(const ImpactParameters & b,
		  const DipoleState & sl, const DipoleState & sr) const {
  nIAccepted = 0;
  nIBelowCut = 0;
  nIPropFail = 0;
  nIKineFail = 0;
  nIOrdering = 0;
  DipoleInteraction::List ret;
  vector<tDipolePtr> dl;
  sl.extract(back_inserter(dl));
  vector<tDipolePtr> dr;
  sr.extract(back_inserter(dr));
  //dont save the lowest fij. Otherwise the FList for LHC PbPb will take too much memory.
  double cutoff = min(0.00000001, 1.0/double(dl.size()*dr.size()));
  if ( dl.size()*dr.size() < 1000000 ) cutoff = 0.0000000001;
  int total = 0;
  int skipped = 0;
  for ( int i = 0, N = dl.size(); i < N; ++i ) {
    for ( int j = 0, M = dr.size(); j < M; ++j ) {
      if ( (M*i + j)%8 == 0 ) UseRandom::rnd();
      /*** TODO: WHY IS VETO FALSE - because we don't want to check
	   kinematics if below cut ***/
      DipoleInteraction di = fij(b, *dl[i], *dr[j], false);
      total++;
      if ( di.f2 > cutoff && !kinematicsVeto(di) ) {
	ret.insert(di);
	++nIAccepted;
      } else {
	skipped++;
	ret.insert(di);
	if ( di.f2 <= cutoff ) ++nIBelowCut;
	else switch ( di.status ) {
	  case DipoleInteraction::PROPFAIL:
	    ++nIPropFail;
	    break;
	  case DipoleInteraction::KINEFAIL:
	    ++nIKineFail;
	    break;
	  case DipoleInteraction::ORDERING:
	    ++nIOrdering;
	    break;
	  case DipoleInteraction::UNKNOWN:
	  case DipoleInteraction::ACCEPTED:
	    break;
	  }
      }
    }
  }
  return ret;
}

DipoleXSec::InteractionRecoil
DipoleXSec::recoil(const DipoleInteraction & di) const {
  InteractionRecoil ret;
  if ( di.norec ) return ret;
  if ( di.ints.first == di.dips.first->partons().first )
    ret.first.first = di.rec;
  else
    ret.first.second = di.rec;
  if ( di.ints.second == di.dips.second->partons().first )
    ret.second.first = di.b->invRotatePT(-di.rec);
  else
    ret.second.second = di.b->invRotatePT(-di.rec);
  return ret;
}

DipoleXSec::InteractionRecoil
DipoleXSec::recoil(const pair<tPartonPtr, tPartonPtr> left,
		   const pair<tPartonPtr, tPartonPtr> right,
		   const ImpactParameters &  b,
		   pair<pair<bool, bool>, pair<bool, bool> > doesInt) const {
  return recoil(left, right, b, doesInt, int0Partons(left.first, left.second,
						     right.first, right.second, b));
}

DipoleXSec::InteractionRecoil
DipoleXSec::recoil(const pair<tPartonPtr, tPartonPtr> left,
		   const pair<tPartonPtr, tPartonPtr> right,
		   const ImpactParameters &  b,
		   pair<pair<bool, bool>, pair<bool, bool> > doesInt,
		   pair<bool, bool> ints) const {
  tPartonPtr p1 = left.first;
  tPartonPtr p2 = left.second;
  tPartonPtr p3 = right.first;
  tPartonPtr p4 = right.second;
  double pTScale = Current<DipoleEventHandler>()->emitter().pTScale();

  if ( theInteraction == 2 ) {  //swing recoil, only new dips
  TransverseMomentum rec14 = -pTScale*b.difference(p1->position(), p4->position())/
    ( b.dist2(*p1,*p4) );
  TransverseMomentum rec23 = -pTScale*b.difference(p2->position(), p3->position())/
    ( b.dist2(*p2,*p3) );

  return make_pair(make_pair(rec14, rec23),
		   make_pair(-rec23, -rec14));
 }

  //4 parton-parton recoils (full diagonals)
  if ( theInteraction == 1 ) {
    TransverseMomentum rec13 = -pTScale*b.difference(p1->position(), p3->position())/
      (b.dist2(*p1,*p3) );
    TransverseMomentum rec14 = -pTScale*b.difference(p1->position(), p4->position())/
      (b.dist2(*p1,*p4) );
    TransverseMomentum rec23 = -pTScale*b.difference(p2->position(), p3->position())/
      (b.dist2(*p2,*p3) );
    TransverseMomentum rec24 = -pTScale*b.difference(p2->position(), p4->position())/
      (b.dist2(*p2,*p4) );
    
    return make_pair(make_pair(rec14 + rec13, rec23 + rec24),
		     make_pair(b.invRotatePT(-rec23 - rec13),
			       b.invRotatePT(-rec14 - rec24)));
  }

  //2 parton-parton recoils (no diagonals)
  if ( theInteraction == 3 ) {
    TransverseMomentum rec13 = -pTScale*b.difference(p1->position(), p3->position())/
      (b.dist2(*p1,*p3) );
    TransverseMomentum rec14 = -pTScale*b.difference(p1->position(), p4->position())/
      (b.dist2(*p1,*p4) );
    TransverseMomentum rec23 = -pTScale*b.difference(p2->position(), p3->position())/
      (b.dist2(*p2,*p3) );
    TransverseMomentum rec24 = -pTScale*b.difference(p2->position(), p4->position())/
      (b.dist2(*p2,*p4) );
    
    if ( doesInt.first.first && doesInt.first.second &&
	 doesInt.second.first && doesInt.second.second )
      return make_pair(make_pair(rec14, rec23),
		       make_pair(b.invRotatePT(-rec23), b.invRotatePT(-rec14)));
    else {
      TransverseMomentum rec11 = TransverseMomentum();
      if ( doesInt.first.first && doesInt.second.first ) rec11 += rec13;
      if ( doesInt.first.first && doesInt.second.second ) rec11 += rec14;
      TransverseMomentum rec12 = TransverseMomentum();
      if ( doesInt.first.second && doesInt.second.first ) rec12 += rec23;
      if ( doesInt.first.second && doesInt.second.second ) rec12 += rec24;
      TransverseMomentum rec21 = TransverseMomentum();
      if ( doesInt.second.first && doesInt.first.first ) rec21 -= rec13;
      if ( doesInt.second.first && doesInt.first.second ) rec21 -= rec23;
      TransverseMomentum rec22 = TransverseMomentum();
      if ( doesInt.second.second && doesInt.first.first ) rec22 -= rec14;
      if ( doesInt.second.second && doesInt.first.second ) rec22 -= rec24;
      return make_pair(make_pair(rec11, rec12),
		       make_pair(b.invRotatePT(rec21), b.invRotatePT(rec22)));
    }
  }

  if ( theInteraction == 0 ) {   //dip-dip recoil

    TransverseMomentum rec13 = -pTScale*b.difference(p1->position(), p3->position())/
      (b.dist2(*p1,*p3) );
    TransverseMomentum rec14 = -pTScale*b.difference(p1->position(), p4->position())/
      (b.dist2(*p1,*p4) );
    TransverseMomentum rec23 = -pTScale*b.difference(p2->position(), p3->position())/
      (b.dist2(*p2,*p3) );
    TransverseMomentum rec24 = -pTScale*b.difference(p2->position(), p4->position())/
      (b.dist2(*p2,*p4) );

    /*** TODO: WHAT IS THIS? ***/
    if ( !Current<DipoleEventHandler>()->eventFiller().compat ) {
      if ( p1->oY() + p3->oY() > p2->oY() + p4->oY() ) rec24 = TransverseMomentum();
      else rec13 = TransverseMomentum();
    }

    // if ( p1->flavour() == ParticleID::g ) {
    //   rec13 /= 2.0;
    //   rec14 /= 2.0;
    // }
    // if ( p2->flavour() == ParticleID::g ) {
    //   rec23 /= 2.0;
    //   rec24 /= 2.0;
    // }
    // if ( p3->flavour() == ParticleID::g ) {
    //   rec13 /= 2.0;
    //   rec23 /= 2.0;
    // }
    // if ( p4->flavour() == ParticleID::g ) {
    //   rec14 /= 2.0;
    //   rec24 /= 2.0;
    // }

    TransverseMomentum rec11 = TransverseMomentum();
    if ( doesInt.first.first && doesInt.second.first ) rec11 += rec13;
    if ( doesInt.first.first && doesInt.second.second ) rec11 += rec14;
    TransverseMomentum rec12 = TransverseMomentum();
    if ( doesInt.first.second && doesInt.second.first ) rec12 += rec23;
    if ( doesInt.first.second && doesInt.second.second ) rec12 += rec24;
    TransverseMomentum rec21 = TransverseMomentum();
    if ( doesInt.second.first && doesInt.first.first ) rec21 -= rec13;
    if ( doesInt.second.first && doesInt.first.second ) rec21 -= rec23;
    TransverseMomentum rec22 = TransverseMomentum();
    if ( doesInt.second.second && doesInt.first.first ) rec22 -= rec14;
    if ( doesInt.second.second && doesInt.first.second ) rec22 -= rec24;

    return make_pair(make_pair(rec11, rec12),
    		     make_pair(b.invRotatePT(rec21), b.invRotatePT(rec22)));

  }

  return make_pair(make_pair(TransverseMomentum(), TransverseMomentum()),
		   make_pair(TransverseMomentum(), TransverseMomentum()));
}

DipoleXSec::RealInteraction
DipoleXSec::initialiseInteraction(const pair<DipolePtr, DipolePtr> inter,
				  RealPartonStatePtr lrs, RealPartonStatePtr rrs,
				  pair<pair<bool, bool>, pair<bool, bool> > doesInt,
				  const ImpactParameters & b) const {
  RealInteraction ret;

  ret.lrs = lrs;
  ret.rrs = rrs;

  ret.d1 = inter.first;
  ret.d2 = inter.second;

  if ( doesInt.first.first ) ret.p11 = lrs->getReal(ret.d1->partons().first);
  else ret.p11 = RealPartonPtr();
  if ( doesInt.first.second ) ret.p12 = lrs->getReal(ret.d1->partons().second);
  else ret.p12 = RealPartonPtr();
  if ( doesInt.second.first ) ret.p21 = rrs->getReal(inter.second->partons().first);
  else ret.p21 = RealPartonPtr();
  if ( doesInt.second.second ) ret.p22 = rrs->getReal(inter.second->partons().second);
  else ret.p22 = RealPartonPtr();

  if ( ret.p11 && ret.p11->fluct != -1 && ret.p12 && ret.p12->fluct == ret.p11->fluct )
    lrs->splitFluct(ret.p11, ret.p12);
  if ( ret.p21 && ret.p21->fluct != -1 && ret.p22 && ret.p22->fluct == ret.p21->fluct )
    rrs->splitFluct(ret.p21, ret.p22);

  ret.range11 = sqrt(min(min(inter.first->partons().first->dist2
			     (*inter.second->partons().first),
			     inter.first->partons().first->dist2
			     (*inter.second->partons().second)),
			 inter.first->partons().first->dist2
			 (*inter.first->partons().second)/4.0));
  ret.range12 = sqrt(min(min(inter.first->partons().second->dist2
			     (*inter.second->partons().first),
			     inter.first->partons().second->dist2
			     (*inter.second->partons().second)),
			 inter.first->partons().second->dist2
			 (*inter.first->partons().first)/4.0));
  ret.range21 = sqrt(min(min(inter.second->partons().first->dist2
			     (*inter.first->partons().first),
			     inter.second->partons().first->dist2
			     (*inter.first->partons().second)),
			 inter.second->partons().first->dist2
			 (*inter.second->partons().second)/4.0));
  ret.range22 = sqrt(min(min(inter.second->partons().second->dist2
			     (*inter.first->partons().first),
			     inter.second->partons().second->dist2
			     (*inter.first->partons().second)),
			 inter.second->partons().second->dist2
			 (*inter.second->partons().first)/4.0));

  if ( theInteraction == 3 ) {
    ret.range11 = sqrt(min(inter.first->partons().first->dist2
			   (*inter.second->partons().second),
			   inter.first->partons().first->dist2
			   (*inter.first->partons().second)/4.0));
    ret.range12 = sqrt(min(inter.first->partons().second->dist2
			   (*inter.second->partons().first),
			   inter.first->partons().second->dist2
			   (*inter.first->partons().first)/4.0));
    ret.range21 = sqrt(min(inter.second->partons().first->dist2
			   (*inter.first->partons().second),
			   inter.second->partons().first->dist2
			   (*inter.second->partons().second)/4.0));
    ret.range22 = sqrt(min(inter.second->partons().second->dist2
			   (*inter.first->partons().first),
			   inter.second->partons().second->dist2
			   (*inter.second->partons().first)/4.0));
  }

  if ( theInteraction == 0 ) {
    pair<bool, bool> ints;
    ints = int0Partons(ret.d1->partons().first, ret.d1->partons().second,
		       inter.second->partons().first, inter.second->partons().second, b);

    InvEnergy2 range2;
    if ( ints.first && ints.second ) range2 = b.dist2(*ret.p11->theParton,*ret.p21->theParton);
    if ( !ints.first && ints.second ) range2 = b.dist2(*ret.p12->theParton,*ret.p21->theParton);
    if ( ints.first && !ints.second ) range2 = b.dist2(*ret.p11->theParton,*ret.p22->theParton);
    if ( !ints.first && !ints.second ) range2 = b.dist2(*ret.p12->theParton,*ret.p22->theParton);
    InvEnergy range = sqrt(range2);

    ret.range11 = min(range, ret.d1->size()/2.0);
    ret.range12 = min(range, ret.d1->size()/2.0);
    ret.range21 = min(range, ret.d2->size()/2.0);
    ret.range22 = min(range, ret.d2->size()/2.0);

    if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 ) {
      ret.range11 = ret.d1->size()/2.0;
      ret.range12 = ret.d1->size()/2.0; 
      ret.range21 = ret.d2->size()/2.0;
      ret.range22 = ret.d2->size()/2.0;
    }
  }

  InvEnergy2 rr11 = (ret.p11 && ret.p21) ? ret.p11->theParton->dist2(*ret.p21->theParton):ZERO;
  InvEnergy2 rr12 = (ret.p11 && ret.p22) ? ret.p11->theParton->dist2(*ret.p22->theParton):ZERO;
  InvEnergy2 rr21 = (ret.p12 && ret.p21) ? ret.p12->theParton->dist2(*ret.p21->theParton):ZERO;
  InvEnergy2 rr22 = (ret.p12 && ret.p22) ? ret.p12->theParton->dist2(*ret.p22->theParton):ZERO;

  InvEnergy2 rm11 = (ret.p11 && ret.p11->mothers.first ?
		    ret.p11->theParton->dist2(*ret.p11->mothers.first->theParton):ZERO);
  InvEnergy2 rm12 = (ret.p12 && ret.p12->mothers.first ?
		    ret.p12->theParton->dist2(*ret.p12->mothers.first->theParton):ZERO);
  InvEnergy2 rm21 = (ret.p21 && ret.p21->mothers.first ?
		    ret.p21->theParton->dist2(*ret.p21->mothers.first->theParton):ZERO);
  InvEnergy2 rm22 = (ret.p22 && ret.p22->mothers.first ?
		    ret.p22->theParton->dist2(*ret.p22->mothers.first->theParton):ZERO);

  ret.max11 = rr11 < rm11 && rr11 < rm21;
  ret.max12 = rr12 < rm11 && rr12 < rm22;
  ret.max21 = rr21 < rm12 && rr21 < rm21;
  ret.max22 = rr22 < rm12 && rr22 < rm22;

  Energy2 den = ((ret.p11 ? 1./rr11:ZERO) + (ret.p12 ? 1./rr12:ZERO)
		 + (ret.p21 ? 1./rr21:ZERO) + (ret.p22 ? 1./rr22:ZERO));
  ret.P11 = ret.p11 ? (1./rr11)/den:ZERO;
  ret.P12 = ret.p12 ? (1./rr12)/den:ZERO;
  ret.P21 = ret.p21 ? (1./rr21)/den:ZERO;
  ret.P22 = ret.p22 ? (1./rr22)/den:ZERO;

  //look only at the new dipole pairs, the cross pairs get zero weight.
  if ( theInteraction == 3 ) {
    den = (1./rr12 + 1./rr21);
    ret.P12 = (1./rr12)/den;
    ret.P21 = (1./rr21)/den;
    ret.P11 = 0.0;
    ret.P22 = 0.0;
  }

  //set only the selected pair to weight 1, the others to 0.
  if ( theInteraction == 0 ) {
    pair<bool, bool> ints;
    ints = int0Partons(inter.first->partons().first, inter.first->partons().second,
		       inter.second->partons().first, inter.second->partons().second, b);

    ret.P11 = 0.0;
    ret.P22 = 0.0;
    ret.P21 = 0.0;
    ret.P12 = 0.0;

    if ( ints.first && ints.second ) ret.P11 = 1.0;
    if ( ints.first && !ints.second ) ret.P12 = 1.0;
    if ( !ints.first && ints.second ) ret.P21 = 1.0;
    if ( !ints.first && !ints.second ) ret.P22 = 1.0;
  }

  return ret;
}

DipoleXSec::RealInteraction
DipoleXSec::initialiseInteraction(const DipoleInteraction & di,
				  RealPartonStatePtr lrs, RealPartonStatePtr rrs) const {
  RealInteraction ret;
  const ImpactParameters & b = *di.b;

  ret.lrs = lrs;
  ret.rrs = rrs;

  ret.d1 = di.dips.first;
  ret.d2 = di.dips.second;

  if ( ret.d1->partons().first == di.ints.first ) ret.p11 = lrs->getReal(di.ints.first);
  if ( ret.d1->partons().second == di.ints.first ) ret.p12 = lrs->getReal(di.ints.first);
  if ( ret.d2->partons().first == di.ints.second ) ret.p21 = rrs->getReal(di.ints.second);
  if ( ret.d2->partons().second == di.ints.second ) ret.p22 = rrs->getReal(di.ints.second);

  if ( ret.p11 && ret.p11->fluct != -1 && ret.p12 && ret.p12->fluct == ret.p11->fluct )
    lrs->splitFluct(ret.p11, ret.p12);
  if ( ret.p21 && ret.p21->fluct != -1 && ret.p22 && ret.p22->fluct == ret.p21->fluct )
    rrs->splitFluct(ret.p21, ret.p22);

  ret.range11 = sqrt(min(min(ret.d1->partons().first->dist2
			     (*ret.d2->partons().first),
			     ret.d1->partons().first->dist2
			     (*ret.d2->partons().second)),
			 ret.d1->partons().first->dist2
			 (*ret.d1->partons().second)/4.0));
  ret.range12 = sqrt(min(min(ret.d1->partons().second->dist2
			     (*ret.d2->partons().first),
			     ret.d1->partons().second->dist2
			     (*ret.d2->partons().second)),
			 ret.d1->partons().second->dist2
			 (*ret.d1->partons().first)/4.0));
  ret.range21 = sqrt(min(min(ret.d2->partons().first->dist2
			     (*ret.d1->partons().first),
			     ret.d2->partons().first->dist2
			     (*ret.d1->partons().second)),
			 ret.d2->partons().first->dist2
			 (*ret.d2->partons().second)/4.0));
  ret.range22 = sqrt(min(min(ret.d2->partons().second->dist2
			     (*ret.d1->partons().first),
			     ret.d2->partons().second->dist2
			     (*ret.d1->partons().second)),
			 ret.d2->partons().second->dist2
			 (*ret.d2->partons().first)/4.0));

  pair<bool, bool> ints(ret.d1->partons().first == di.ints.first,
			ret.d2->partons().first == di.ints.second);
  
  InvEnergy2 range2 = b.dist2(*di.ints.first, *di.ints.second);
  InvEnergy range = sqrt(range2);

  ret.range11 = min(range, ret.d1->size()/2.0);
  ret.range12 = min(range, ret.d1->size()/2.0);
  ret.range21 = min(range, ret.d2->size()/2.0);
  ret.range22 = min(range, ret.d2->size()/2.0);

  if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 ) {
    ret.range11 = ret.d1->size()/2.0;
    ret.range12 = ret.d1->size()/2.0; 
    ret.range21 = ret.d2->size()/2.0;
    ret.range22 = ret.d2->size()/2.0;
  }

  InvEnergy2 rr11 = (ret.p11 && ret.p21) ? ret.p11->theParton->dist2(*ret.p21->theParton):ZERO;
  InvEnergy2 rr12 = (ret.p11 && ret.p22) ? ret.p11->theParton->dist2(*ret.p22->theParton):ZERO;
  InvEnergy2 rr21 = (ret.p12 && ret.p21) ? ret.p12->theParton->dist2(*ret.p21->theParton):ZERO;
  InvEnergy2 rr22 = (ret.p12 && ret.p22) ? ret.p12->theParton->dist2(*ret.p22->theParton):ZERO;

  InvEnergy2 rm11 = (ret.p11 && ret.p11->mothers.first ?
		    ret.p11->theParton->dist2(*ret.p11->mothers.first->theParton):ZERO);
  InvEnergy2 rm12 = (ret.p12 && ret.p12->mothers.first ?
		    ret.p12->theParton->dist2(*ret.p12->mothers.first->theParton):ZERO);
  InvEnergy2 rm21 = (ret.p21 && ret.p21->mothers.first ?
		    ret.p21->theParton->dist2(*ret.p21->mothers.first->theParton):ZERO);
  InvEnergy2 rm22 = (ret.p22 && ret.p22->mothers.first ?
		    ret.p22->theParton->dist2(*ret.p22->mothers.first->theParton):ZERO);

  ret.max11 = rr11 < rm11 && rr11 < rm21;
  ret.max12 = rr12 < rm11 && rr12 < rm22;
  ret.max21 = rr21 < rm12 && rr21 < rm21;
  ret.max22 = rr22 < rm12 && rr22 < rm22;

  Energy2 den = ((ret.p11 ? 1./rr11:ZERO) + (ret.p12 ? 1./rr12:ZERO)
		 + (ret.p21 ? 1./rr21:ZERO) + (ret.p22 ? 1./rr22:ZERO));
  ret.P11 = ret.p11 ? (1./rr11)/den:ZERO;
  ret.P12 = ret.p12 ? (1./rr12)/den:ZERO;
  ret.P21 = ret.p21 ? (1./rr21)/den:ZERO;
  ret.P22 = ret.p22 ? (1./rr22)/den:ZERO;

  //look only at the new dipole pairs, the cross pairs get zero weight.
  if ( theInteraction == 3 ) {
    den = (1./rr12 + 1./rr21);
    ret.P12 = (1./rr12)/den;
    ret.P21 = (1./rr21)/den;
    ret.P11 = 0.0;
    ret.P22 = 0.0;
  }

  //set only the selected pair to weight 1, the others to 0.

  ret.P11 = 0.0;
  ret.P22 = 0.0;
  ret.P21 = 0.0;
  ret.P12 = 0.0;

  if ( ints.first && ints.second ) ret.P11 = 1.0;
  if ( ints.first && !ints.second ) ret.P12 = 1.0;
  if ( !ints.first && ints.second ) ret.P21 = 1.0;
  if ( !ints.first && !ints.second ) ret.P22 = 1.0;

  return ret;
}

pair<pair<bool, bool>, pair<bool, bool> >
DipoleXSec::doesInt(const pair<tPartonPtr, tPartonPtr> left,
		    const pair<tPartonPtr, tPartonPtr> right,
		    const ImpactParameters & b) const {
  //decide which of the 4 partons will end up on shell

  // if ( Current<DipoleEventHandler>()->eventFiller().compat ) {
  //   pair<bool,bool> int0 =
  //     int0Partons(left.first, left.second, right.first, right.second, b);
  //   return make_pair(make_pair(int0.first, !int0.first),
  // 		     make_pair(int0.second, !int0.second));
  // } else {
    return make_pair(make_pair(true, true), make_pair(true, true));
  // }

  if ( Current<DipoleEventHandler>()->eventFiller().singleMother() ) {
    //when single mother not all are chosen
    pair<pair<bool, bool>, pair<bool, bool>  > ret =
      make_pair(make_pair(false, false), make_pair(false, false));


    if ( theInteraction == 0 ) {
      //if interaction 0, then always keep the same partons used in fij
      pair<bool, bool> ints;
      ints = int0Partons(left.first, left.second, right.first, right.second, b);

      if ( ints.first ) ret.first.first = true;
      else ret.first.second = true;
      if ( ints.second ) ret.second.first =true;
      else ret.second.second = true;
    }
    else {
      //select which two partons the gluons are exchanged between,
      //weighted by 1/distance^2
      InvEnergy2 r11 = b.dist2(*left.first, *right.first);
      InvEnergy2 r12 = b.dist2(*left.first, *right.second);
      InvEnergy2 r21 = b.dist2(*left.second, *right.first);
      InvEnergy2 r22 = b.dist2(*left.second, *right.second);

      Selector<int, double> sel;
      sel.insert(1./r11/(1./r11 + 1./r12 + 1./r21 + 1./r22), 1);
      sel.insert(1./r12/(1./r11 + 1./r12 + 1./r21 + 1./r22), 2);
      sel.insert(1./r21/(1./r11 + 1./r12 + 1./r21 + 1./r22), 3);
      sel.insert(1./r22/(1./r11 + 1./r12 + 1./r21 + 1./r22), 4);

      double dummy = (r11 + r22 + r12 + r21)*sqr(GeV);
      double pseudoRnd = 10000.0*dummy - int(10000.0*dummy);

      switch ( sel[pseudoRnd] ) {
      case 1:
	ret.first.first = true;
	ret.second.first = true;
	break;
      case 2:
	ret.first.first = true;
	ret.second.second = true;
	break;
      case 3:
	ret.first.second = true;
	ret.second.first = true;
	break;
      case 4:
	ret.first.second = true;
	ret.second.second = true;
	break;
      }
    }

    //check for swinged emission, in which case both partons are real
    if ( !(left.first->parents().second == left.second ||
	   left.second->parents().first == left.first) ) {
      ret.first.first = true;
      ret.first.second = true;
    }
    if ( !(right.first->parents().second == right.second ||
	   right.second->parents().first == right.first) ) {
      ret.second.first = true;
      ret.second.second = true;
    }
    return ret;
  }

  //if not single mother, all partons are on shell
  return make_pair(make_pair(true, true), make_pair(true, true));
}

pair<bool, bool> DipoleXSec::int0Partons(tcPartonPtr p11, tcPartonPtr p12,
					 tcPartonPtr p21, tcPartonPtr p22, 
					 const ImpactParameters & b) const {
  //chose partons (pseudo-)randomly
  InvEnergy2 rr11 = b.dist2(*p11, *p21);
  InvEnergy2 rr22 = b.dist2(*p12, *p22);
  double dummy = (rr11 + rr22)*GeV2;
  while ( dummy > 100 ) dummy /= 10;

  if ( Current<DipoleEventHandler>()->fudgeME() > 1 ) {
    pair<bool, bool> ret = make_pair(true, true);
    if ( p11->dipoles().second->colour() == p21->dipoles().second->colour() )
      ret.second = false;
    if ( 1000000.0*dummy - long(1000000.0*dummy) > 0.5 ) {
      ret.first = !ret.first;
      ret.second = !ret.second;
    }
    return ret;
  }

  bool firstDipole  = (  1000000.0*dummy -  long(1000000.0*dummy) > 0.5 );
  bool secondDipole = ( 10000000.0*dummy - long(10000000.0*dummy) > 0.5 );

  return make_pair(firstDipole, secondDipole);
}

void DipoleXSec::updateMomenta(RealInteraction * i) const {
  if ( i->p11 ) {
    pair<Energy, Energy> pm11 = i->p11->effectivePlusMinus(i->range11, true);
    i->effectivePlus11 = pm11.first;
    i->effectiveMinus11 = pm11.second;
  }
  if ( i->p12 ) {
    pair<Energy, Energy> pm12 = i->p12->effectivePlusMinus(i->range12, false);
    i->effectivePlus12 = pm12.first;
    i->effectiveMinus12 = pm12.second;
  }
  if ( i->p21 ) {
    pair<Energy, Energy> pm21 = i->p21->effectivePlusMinus(i->range21, true);
    i->effectivePlus21 = pm21.first;
    i->effectiveMinus21 = pm21.second;
  }
  if ( i->p22 ) {
    pair<Energy, Energy> pm22 = i->p22->effectivePlusMinus(i->range22, false);
    i->effectivePlus22 = pm22.first;
    i->effectiveMinus22 = pm22.second;
  }
}

void DipoleXSec::doTransverseRecoils(RealInteraction i, InteractionRecoil recs) const {
  if ( Current<DipoleEventHandler>()->eventFiller().mode() != 1 ) {
    i.lrs->totalRecoil += recs.first.first + recs.first.second;
    i.rrs->totalRecoil += recs.second.first + recs.second.second;
    if ( i.p11 ) i.p11->intRecoil += recs.first.first;
    if ( i.p12 ) i.p12->intRecoil += recs.first.second;
    if ( i.p21 ) i.p21->intRecoil += recs.second.first;
    if ( i.p22 ) i.p22->intRecoil += recs.second.second;
  }

  if ( i.p11 ) {
    if ( i.p22 ) i.p22->doEffectiveRecoil(i.p11, i.range11, true, ZERO,  -recs.first.first);
    else i.p21->doEffectiveRecoil(i.p11, i.range11, true, ZERO,  -recs.first.first);
  }
  if ( i.p12 ) {
    if ( i.p22 ) i.p22->doEffectiveRecoil(i.p12, i.range12, false, ZERO,  -recs.first.second);
    else i.p21->doEffectiveRecoil(i.p12, i.range12, false, ZERO,  -recs.first.second);
  }
  if ( i.p21 ) {
    if ( i.p12 ) i.p12->doEffectiveRecoil(i.p21, i.range21, true, ZERO,  -recs.second.first);
    else i.p11->doEffectiveRecoil(i.p21, i.range21, true, ZERO,  -recs.second.first);
  }
  if ( i.p22 ) {
    if ( i.p11 ) i.p11->doEffectiveRecoil(i.p22, i.range22, false, ZERO,  -recs.second.second);
    else i.p12->doEffectiveRecoil(i.p22, i.range22, false, ZERO,  -recs.second.second);
  }

}

pair<double, double> DipoleXSec::findBoosts(RealInteraction i) const {
  Energy neededMinus = ZERO; 
  if ( checkOffShell ) neededMinus += i.lrs->neededValenceMinus();
  if ( i.p11 ) neededMinus += i.p11->effectiveGiveMinus(i.range11, true);
  if ( i.p12 ) neededMinus += i.p12->effectiveGiveMinus(i.range12, false);
  Energy neededPlus = ZERO;
  if ( checkOffShell ) neededPlus += i.rrs->neededValenceMinus();
  if ( i.p21 ) neededPlus += i.p21->effectiveGiveMinus(i.range21, true);
  if ( i.p22 ) neededPlus += i.p22->effectiveGiveMinus(i.range22, false);

  Energy intPlus1 = ZERO;
  if ( i.p11 ) intPlus1 += i.effectivePlus11;
  if ( i.p12 ) intPlus1 += i.effectivePlus12;
  Energy intPlus2 = ZERO;
  if ( i.p21 ) intPlus2 += i.p21->inRangeMinus(i.range21, true);
  if ( i.p22 ) intPlus2 += i.p22->inRangeMinus(i.range22, false);

  Energy intMinus1 = ZERO;
  if ( i.p11 ) intMinus1 += i.p11->inRangeMinus(i.range11, true);
  if ( i.p12 ) intMinus1 += i.p12->inRangeMinus(i.range12, false);
  Energy intMinus2 = ZERO;
  if ( i.p21 ) intMinus2 += i.effectivePlus21;
  if ( i.p22 ) intMinus2 += i.effectivePlus22;

  Energy evoPlus2 = neededPlus - intPlus2;
  Energy evoMinus1 = neededMinus - intMinus1;

  pair<double, double> boosts = findBoosts(intPlus1, intPlus2, intMinus1, intMinus2, evoPlus2, evoMinus1);

  if ( boosts.first == 0.0 || isnan(boosts.first) || isnan(boosts.second) ) {
    return make_pair(0.0,0.0);
  }
  if ( boosts.second < 0.0 || boosts.first < 0.0 ) {
    //can be over 1.0 if the valenceminus is more than needed
    return make_pair(0.0,0.0);;
  }

  return boosts;
}

void DipoleXSec::doBoosts(RealInteraction i, pair<double, double> boosts) const {
  doBoost(i.p11, i.range11, i.p12, i.range12, boosts.first);
  doBoost(i.p21, i.range21, i.p22, i.range22, boosts.second);
}

void DipoleXSec::reduceRecoil(RealInteraction RI, InteractionRecoil recs) const {
  //check if the interacting partons have passed each other
  while ( max((RI.p11 ? RI.p11->y:RI.p12->y), (RI.p12 ? RI.p12->y:RI.p11->y)) >
	  -max((RI.p21 ? RI.p21->y:RI.p22->y), (RI.p22 ? RI.p22->y:RI.p21->y)) ) {
    //and if the recoil is larger than 1 GeV
    if ( max(max(recs.first.first.pt(), recs.first.second.pt()),
	     max(recs.second.first.pt(), recs.second.second.pt())) < 1*GeV ) {
      break;
    }
    //then undo the recoil and redo 90% of it.
    recs = make_pair(make_pair(-recs.first.first, -recs.first.second),
		     make_pair(-recs.second.first, -recs.second.second));
    doTransverseRecoils(RI, recs);
    recs = make_pair(make_pair(-recs.first.first*0.9, -recs.first.second*0.9),
		     make_pair(-recs.second.first*0.9, -recs.second.second*0.9));
    doTransverseRecoils(RI, recs);
  }
}

bool DipoleXSec::doInteraction(InteractionRecoil recs, const FList::const_iterator inter,
			       RealPartonStatePtr lrs, RealPartonStatePtr rrs,
			       pair<pair<bool, bool>, pair<bool, bool> > doesInt,
			       const ImpactParameters & b) const {
  static DebugItem checkkinematics("DIPSY::CheckKinematics", 6);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found at start of interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found at start of interaction!" << Exception::warning;

  RealInteraction RI = initialiseInteraction(inter->second, lrs, rrs, doesInt, b);

  doTransverseRecoils(RI, recs);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after transverse recoil in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after transverse recoil in interaction!" << Exception::warning;

  if ( theRecoilReduction ) {
    reduceRecoil(RI, recs);
  }

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after recoil reduction in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after recoil reduction in interaction!" << Exception::warning;

  updateMomenta(& RI);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after update in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after update in interaction!" << Exception::warning;

  pair<double, double> boosts = findBoosts(RI);
  if ( boosts.first == 0.0 ) {
    return false;
  }

  doBoosts(RI, boosts);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after boost in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after boost in interaction!" << Exception::warning;

  if ( ordered(RI, recs, b) ) {
    setOnShell(RI);
    if ( checkkinematics && lrs->checkForNegatives() )
      Throw<InteractionKinematicException>()
	<< "negatives found at end of interaction!" << Exception::warning;
    if ( checkkinematics && rrs->checkForNegatives() )
      Throw<InteractionKinematicException>()
	<< "negatives found at end of interaction!" << Exception::warning;
    return true;
  }
  return false;
}

bool DipoleXSec::doInteraction(InteractionRecoil recs,
			       const DipoleInteraction::List::const_iterator inter,
			       RealPartonStatePtr lrs, RealPartonStatePtr rrs,
			       pair<pair<bool, bool>, pair<bool, bool> > doesInt) const {
  const ImpactParameters & b = *(inter->b);

  static DebugItem checkkinematics("DIPSY::CheckKinematics", 6);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found at start of interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found at start of interaction!" << Exception::warning;

  RealInteraction RI = initialiseInteraction(*inter, lrs, rrs);

  doTransverseRecoils(RI, recs);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after transverse recoil in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after transverse recoil in interaction!" << Exception::warning;

  if ( theRecoilReduction ) {
    reduceRecoil(RI, recs);
  }

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after recoil reduction in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after recoil reduction in interaction!" << Exception::warning;

  updateMomenta(& RI);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after update in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after update in interaction!" << Exception::warning;

  pair<double, double> boosts = findBoosts(RI);
  if ( boosts.first == 0.0 ) {
    return false;
  }

  doBoosts(RI, boosts);

  if ( checkkinematics && lrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after boost in interaction!" << Exception::warning;
  if ( checkkinematics && rrs->checkForNegatives() )
    Throw<InteractionKinematicException>()
      << "negatives found after boost in interaction!" << Exception::warning;

  if ( ordered(RI, recs, b) ) {
    setOnShell(RI);
    if ( checkkinematics && lrs->checkForNegatives() )
      Throw<InteractionKinematicException>()
	<< "negatives found at end of interaction!" << Exception::warning;
    if ( checkkinematics && rrs->checkForNegatives() )
      Throw<InteractionKinematicException>()
	<< "negatives found at end of interaction!" << Exception::warning;
    return true;
  }
  return false;

}


bool DipoleXSec::ordered(RealInteraction i, InteractionRecoil recs,
			 const ImpactParameters & b) const {
  //this options means only momentum conservation, and that is already checked by
  //finding a boost. no extra ordering wanted.
  if ( theIntOrdering == 2 ) {
    return true;
  }

  double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();

  bool ordered = true;
  bool both = Current<DipoleEventHandler>()->emitter().bothOrderedFS();
  //if not a local pt max, check p+- ordering
  if ( theIntOrdering == 0 ) {
    if ( theInteraction == 0 ) {
      //check ordering only for the key partons
      pair<bool, bool> ints0 = int0Partons(i.d1->partons().first, i.d1->partons().second,
					   i.d2->partons().first, i.d2->partons().second, b);
      Energy effMinus1 = (ints0.first ? i.effectiveMinus11:i.effectiveMinus12);
      Energy effMinus2 = (ints0.second ? i.effectiveMinus21:i.effectiveMinus22);
      Energy effPlus1 = (ints0.first ? i.effectivePlus11:i.effectivePlus12);
      Energy effPlus2 = (ints0.second ? i.effectivePlus21:i.effectivePlus22);

      if ( effPlus1*PSInf < effMinus2 )  ordered = false;
      if ( effPlus2*PSInf < effMinus1 )  ordered = false;
    }
    else {
      if ( !i.max11 && i.p11 && i.p21 ) {
	if ( i.effectivePlus11*PSInf < i.effectiveMinus21*(both ? 1.0:i.P11) )  ordered = false;
	if ( (both ? 1.0:i.P11)*i.effectiveMinus11/PSInf > i.effectivePlus21 )  ordered = false;
      }
      if ( !i.max12 && i.p11 && i.p22 ) {
	if ( i.effectivePlus11*PSInf < i.effectiveMinus22*(both ? 1.0:i.P12) )  ordered = false;
	if ( (both ? 1.0:i.P12)*i.effectiveMinus11/PSInf > i.effectivePlus22 )  ordered = false;
      }
      if ( !i.max21 && i.p12 && i.p21 ) {
	if ( i.effectivePlus12*PSInf < i.effectiveMinus21*(both ? 1.0:i.P21) )  ordered = false;
	if ( (both ? 1.0:i.P21)*i.effectiveMinus12/PSInf > i.effectivePlus21 )  ordered = false;
      }
      if ( !i.max22 && i.p12 && i.p22 ) {
	if ( i.effectivePlus12*PSInf < i.effectiveMinus22*(both ? 1.0:i.P22) )  ordered = false;
	if ( (both ? 1.0:i.P22)*i.effectiveMinus12/PSInf > i.effectivePlus22 )  ordered = false;
      }
    }
  }

  if ( theIntOrdering == 1 && theInteraction == 0 ) {
    double PMOrd = Current<DipoleEventHandler>()->emitter().PMinusOrdering();

    //check ordering only for the key partons
    pair<bool, bool> ints0 = int0Partons(i.d1->partons().first, i.d1->partons().second,
					 i.d2->partons().first, i.d2->partons().second, b);
    Energy effMinus1 = (ints0.first ? i.effectiveMinus11:i.effectiveMinus12);
    Energy effMinus2 = (ints0.second ? i.effectiveMinus21:i.effectiveMinus22);
    Energy rec1 = (ints0.first ? recs.first.first.pt():recs.first.second.pt());
    Energy rec2 = (ints0.second ? recs.second.first.pt():recs.second.second.pt());

    if ( sqr(max(rec1, rec2)*PSInf) < effMinus1*effMinus2*PMOrd ) ordered = false;
  }
  
  if ( (i.p11 && i.p11->searchNegative(i.range11, true)) ||
       (i.p12 && i.p12->searchNegative(i.range12, false)) ||
       (i.p21 && i.p21->searchNegative(i.range21, true)) ||
       (i.p22 && i.p22->searchNegative(i.range22, false)) )
      ordered = false;

  return ordered;
}

void DipoleXSec::setOnShell(RealInteraction i) const {
  if ( i.p11 ) i.p11->effectiveGiveMinus(i.range11, true);
  if ( i.p12 ) i.p12->effectiveGiveMinus(i.range12, false);
  if ( i.p21 ) i.p21->effectiveGiveMinus(i.range21, true);
  if ( i.p22 ) i.p22->effectiveGiveMinus(i.range22, false);
}


bool DipoleXSec::reconnect(tDipolePtr d1, tDipolePtr d2) const {

  if ( d1->children().second && !d1->children().first ) {
    return reconnect(d1->children().second, d2);
  }

  if ( d2->children().second && !d2->children().first ) {
    return reconnect(d1, d2->children().second);
  }
  if ( !d1->children().first && !d2->children().first ) { //none has rescattered
    d1->swingDipole(d2);
    Current<DipoleEventHandler>()->swinger().recombine(*d1);
    interact(*d1->children().first, *d2->children().first);
    //    d1->children().first->interact(*d2->children().first);
    //    d2->children().first->interact(*d1->children().first);
    return true;
  }
  tPartonPtr p11 = d1->partons().first;
  tPartonPtr p12 = d1->partons().second;
  tPartonPtr p21 = d2->partons().first;
  tPartonPtr p22 = d2->partons().second;
  tDipolePtr d11 = p11->dipoles().second;
  tDipolePtr d12 = p12->dipoles().first;
  tDipolePtr d21 = p21->dipoles().second;
  tDipolePtr d22 = p22->dipoles().first;
  tDipolePtr swing1;
  tDipolePtr swing2;
  if ( d11 == d12 ) d1 = d11; //dipole has swinged back
  if ( d21 == d22 ) d2 = d21;

  if ( !d1->children().first && !d2->children().first ) { //both are original dipole
    swing1 = d1;
    swing2 = d2;
  }
  else if ( d1->children().first && !d2->children().first ) { //one dip d1 is rescattering
    tDipolePtr temp = d1;
    d1 = d2;
    d2 = temp;
    p11 = d1->partons().first;
    p12 = d1->partons().second;
    p21 = d2->partons().first;
    p22 = d2->partons().second;
    d11 = p11->dipoles().second;
    d12 = p12->dipoles().first;
    d21 = p21->dipoles().second;
    d22 = p22->dipoles().first;
  }
  if ( !d1->children().first && d2->children().first ) { //one dip d2 is rescattering
    if ( p11 != d21->partons().second && p12 != d22->partons().first ) {
      //no connections, swing with one of the two randomly
      swing1 = d1;
      if ( UseRandom::rnd() > 0.5 ) swing2 = d21;
      else swing2 = d22;
    }
    else if ( p11 == d21->partons().second && p12 != d22->partons().first ) {
      //share one swinged dip, swing with the other
      swing1 = d1;
      swing2 = d22;
    }
    else if ( p11 != d21->partons().second && p12 == d22->partons().first ) {
      //share other swinged dip, swing with the first
      swing1 = d1;
      swing2 = d21;
    }
    else if ( p11 == d21->partons().second && p12 == d22->partons().first ) {
      //share both swinged dip, swing back to original
      swing1 = d21;
      swing2 = d22;
    }
  }

  else if ( d1->children().first && d2->children().first ) { //both dips are rescattering
    bool found = false;
    int i = 0;
    while ( !found && i < 1000 ) {
      if ( UseRandom::rnd() > 0.5 ) swing1 = d11;
      else swing1 = d12;
      if ( UseRandom::rnd() > 0.5 ) swing2 = d21;
      else swing2 = d22;
      if ( swing1 != swing2 && swing1->partons().first != swing2->partons().second &&
	   swing2->partons().first != swing1->partons().second )
	found = true; 
   }
  }
  if ( !swing1 || !swing2 ) {
    d1->dipoleState().diagnosis(true);
    return false;
  }
  swing1->swingDipole(swing2);
  Current<DipoleEventHandler>()->swinger().recombine(*swing1);
  interact(*swing1->children().first, *swing2->children().first);
  // swing1->children().first->interact(*swing2->children().first);
  // swing2->children().first->interact(*swing1->children().first);
  return true;
}

void DipoleXSec::interact(Dipole & d1, Dipole & d2) const {
  d1.interact(d2, partonicInteraction());
  d2.interact(d1, partonicInteraction());
}

vector<pair<DipolePtr, DipolePtr> >
DipoleXSec::getColourExchanges(tRealPartonStatePtr lrs, tRealPartonStatePtr rrs) const {
  vector<pair<DipolePtr, DipolePtr> > ret;
  while ( ret.empty() ) {
    list<tDipolePtr>::iterator rightDip = rrs->interactions.begin();
    for ( list<tDipolePtr>::iterator leftDip = lrs->interactions.begin();
	  leftDip != lrs->interactions.end(); leftDip++, rightDip++ ) {
      if ( unitarize(fij((*leftDip)->partons(), (*rightDip)->partons(),
			 ImpactParameters(), false))*0.99
	   < UseRandom::rnd() || true )
	ret.push_back(make_pair(*leftDip, *rightDip));
    }
  }
  return ret;
}

pair< double, double> DipoleXSec::findBoosts(Energy intPlus1, Energy intPlus2,
					     Energy intMinus1, Energy intMinus2,
					     Energy evoPlus2, Energy evoMinus1) const {
  Energy2 A = (intPlus1 - evoPlus2)*intMinus2;
  Energy2 B = intPlus1*(evoMinus1 + intMinus1 - intMinus2) - 
    evoPlus2*(evoMinus1 - intMinus2) - intPlus2*intMinus2;
  Energy2 C = -intPlus2*(evoMinus1 - intMinus2);
  if ( sqr(B/(2*A)) - C/A < 0.0 ) return make_pair(0.0, 0.0);
  double y = -B/(2*A) + sqrt(sqr(B/(2*A)) - C/A); //is the + sqrt() always the right solution?
  double x = 1.0 - intPlus2/(y*intPlus1) - evoPlus2/intPlus1;
  return make_pair(x, y);
}

void DipoleXSec::doBoost(tRealPartonPtr p1, InvEnergy range1,
			 tRealPartonPtr p2, InvEnergy range2, double x) const {
  //we here have to use the plusweighted recoil to fit the x, y above.
  //Other weights get a lot more complicated equations for the boosts.
  if ( p1 ) {
    pair<Energy, Energy> pm1 = p1->effectivePlusMinus(range1, true);
    p1->doPlusWeightedRecoil(p1, range1, true, (1.0 - x)*pm1.first, TransverseMomentum());
  }
  if ( p2 ) {
    pair<Energy, Energy> pm2 = p2->effectivePlusMinus(range2, false);
    p2->doPlusWeightedRecoil(p2, range2, false, (1.0 - x)*pm2.first, TransverseMomentum());
  }
}

double DipoleXSec::unitarize(double f) const {
  return Math::exp1m(-f);
}

void DipoleXSec::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRMax, InvGeV) << theInteraction << sinFunction << usePartonicInteraction
     << theIntOrdering << theRecoilReduction << checkOffShell;
}

void DipoleXSec::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRMax, InvGeV) >> theInteraction >> sinFunction >> usePartonicInteraction
     >> theIntOrdering >> theRecoilReduction >> checkOffShell;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<DipoleXSec,HandlerBase>
  describeDIPSYDipoleXSec("DIPSY::DipoleXSec", "libAriadne5.so libDIPSY.so");

void DipoleXSec::Init() {

  static ClassDocumentation<DipoleXSec> documentation
    ("There is no documentation for the DipoleXSec class");

  static Parameter<DipoleXSec,InvEnergy> interfaceRMax
    ("RMax",
     "The confinement scale (in iverse GeV). If set to zero, "
     "the value of <interface>DipoleEventHandler::RMax</interface> of the "
     "controlling event handler will be used.",
     &DipoleXSec::theRMax, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<DipoleXSec, int> interfaceInteraction
    ("Interaction",
     "Which interaction to be used. Determines f_{ij} and recoils. "
     "0 is the sinus interaction with 4 identical recoils"
     "1 is the sinus interaction with recoils between all 4 parton pairs"
     "2 is the swing interaction"
     "3 is the sinus interaction with recoils between the 2 parton pairs"
     "                           that get the new dipoles",
     &DipoleXSec::theInteraction, 1, 1, 0, 0,
     true, false, Interface::lowerlim);


  static Switch<DipoleXSec,int> interfaceSinFunction
    ("SinFunction",
     "Determine what approximation to the sine functions in the interaction strength "
     "to use.",
     &DipoleXSec::sinFunction, 0, true, false);
  static SwitchOption interfaceSinFunctionExact
    (interfaceSinFunction,
     "Exact",
     "The actual function.",
     0);
  static SwitchOption interfaceSinFunctionAverage
    (interfaceSinFunction,
     "Average",
     "Average over angle of dipole and of the resulting bessel fuction.",
     1);

  static Switch<DipoleXSec,int> interfaceIntOrdering
    ("IntOrdering",
     "How the real state is found from the virtual cascade. "
     "Speed versus consistency.",
     &DipoleXSec::theIntOrdering, 0, true, false);
  static SwitchOption interfaceIntOrderingDefault
    (interfaceIntOrdering,
     "Default",
     "plus and minus reuired to be ordered after all recoils.",
     0);
  static SwitchOption interfaceIntOrderingAsEvo
    (interfaceIntOrdering,
     "AsEvo",
     "Try to emulate the ordering in the evolution by requesting ordering"
     " on the colliding particle if it would've been part of the cascade."
     " That is, ignore recoils from the colliding cascade.",
     1);
  static SwitchOption interfaceIntOrderingVeryOpen
    (interfaceIntOrdering,
     "VeryOpen",
     "Only checks that there is enough energy to put the interaction recoil "
     "on shell. Does not care about ordering, or setting evolution pt on "
     "shell. Same as Emils code.",
     2);
  static SwitchOption interfaceIntOrderingShadowOpen
    (interfaceIntOrdering,
     "ShadowOpen",
     "Only checks that there is enough energy to put the interaction recoil "
     "on shell.",
     3);
  static SwitchOption interfaceIntOrderingShadowColour
    (interfaceIntOrdering,
     "ShadowColour",
     "Check that all partons on the same colour line are ordered.",
     4);

  static Switch<DipoleXSec,int> interfaceRecoilReduction
    ("RecoilReduction",
     "What to do with large recoils",
     &DipoleXSec::theRecoilReduction, 0, true, false);
  static SwitchOption interfaceRecoilReductionOff
    (interfaceRecoilReduction,
     "Off",
     "Do nothing special.",
     0);
  static SwitchOption interfaceRecoilReductionRapidityOrdered
    (interfaceRecoilReduction,
     "RapidityOrdered",
     "Reduce the recoil until all interacting partons are rapidity ordered with each other.",
     1);


  static Switch<DipoleXSec,bool> interfaceCheckOffShell
    ("CheckOffShell",
     "Make sure there is energy available to put incoming particles on-shell. Only necessary for virtual photons.",
     &DipoleXSec::checkOffShell, true, true, false);
  static SwitchOption interfaceCheckOffShellYes
    (interfaceCheckOffShell,
     "Yes",
     "Do the check",
     true);
  static SwitchOption interfaceCheckOffShellNo
    (interfaceCheckOffShell,
     "No",
     "No check is made.",
     false);


  static Switch<DipoleXSec,bool> interfacePartonicInteraction
    ("PartonicInteraction",
     "Flag determining if only one parton in each dipole is considered interacting or both.",
     &DipoleXSec::usePartonicInteraction, false, true, false);
  static SwitchOption interfacePartonicInteractionDipole
    (interfacePartonicInteraction,
     "Dipole",
     "The both partons in a dipole interact.",
     false);
  static SwitchOption interfacePartonicInteractionParton
    (interfacePartonicInteraction,
     "Parton",
     "Only one parton in a dipole interacts.",
     true);

}

vector<int> DipoleInteraction::ofail(18, 0);
vector<int> DipoleInteraction::o1fail(18, 0);
