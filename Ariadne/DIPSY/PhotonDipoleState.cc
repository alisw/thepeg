// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhotonDipoleState class.
//

#include "PhotonDipoleState.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Current.h"
#include "DipoleEventHandler.h"


#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PhotonDipoleState.tcc"
#endif


using namespace DIPSY;

PhotonDipoleState::
PhotonDipoleState(const DipoleEventHandler & eh, Energy plus, Energy minus,
		  Ptr<PhotonWFInfo>::pointer wfi, double wgt)
  : DipoleState(eh, wfi) {
  PPtr inc = wfi->wf().particle()->produceParticle(lightCone(plus, minus));
  inc->getInfo().push_back(wfi);
  vector<PartonPtr> & partons = theIncoming[inc];
  thePlus = plus;
  theMinus = minus;
  weight(wgt);
  PartonPtr q = new_ptr(Parton());
  partons.push_back(q);
  PartonPtr qb = new_ptr(Parton());
  partons.push_back(qb);
  DipolePtr d = createDipole();
  d->partons(make_pair(q, qb));
  q->dipoles(make_pair(DipolePtr(), d));
  qb->dipoles(make_pair(d, DipolePtr()));
  generateColourIndex(d);

  double phi = UseRandom::rnd()*2.0*Constants::pi;
  Parton::Point point(wfi->r()*cos(phi)/2.0, wfi->r()*sin(phi)/2.0);
  q->position(point);
  qb->position(-point);

  TransverseMomentum pt = q->recoil(qb);
  q->pT(pt);
  q->plus(wfi->z()*plus);
  q->updateYMinus();
  q->oY(q->y());
  q->flavour(wfi->flav());

  qb->pT(-pt);
  qb->plus((1.0 - wfi->z())*plus);
  qb->updateYMinus();
  qb->oY(qb->y());
  qb->flavour(-wfi->flav());

  q->valencePT(q->pT());
  qb->valencePT(qb->pT());
  q->valencePlus(q->plus());
  qb->valencePlus(qb->plus());
  theInitialDipoles.push_back(d);

  //trying to conserve 4-momentum


  theMinusDeficit = - minus + q->minus() + qb->minus();
  // cout << "new state with plus, minus: " << plus/GeV << ", " << minus/GeV << endl;
  // cout << "minusDeficit: " << theMinusDeficit/GeV << endl;


}

PhotonDipoleState::~PhotonDipoleState() {}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeNoPIOClass<PhotonDipoleState,DIPSY::DipoleState>
  describeDIPSYPhotonDipoleState("DIPSY::PhotonDipoleState", "libAriadne5.so libDIPSY.so");


void PhotonDipoleState::Init() {}

