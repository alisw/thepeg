// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleNucleusState class.
//

#include "SimpleNucleusState.h"
#include "SimpleNucleus.h"
#include "DipoleEventHandler.h"
#include "Dipole.h"
#include "Parton.h"
#include "ThePEG/Repository/UseRandom.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SimpleNucleusState.tcc"
#endif


using namespace DIPSY;

SimpleNucleusState::
SimpleNucleusState(const DipoleEventHandler & eh, SimpleNucleus & WF,
		   Energy plus, Energy minus, const vector<Point> & positions,
		   WFInfoPtr wfi)
  : DipoleState(eh, wfi) {

  int isys = 0;
  int np = 0;
  int nn = 0;
  if ( WF.proton() && WF.neutron() ) {
    np = abs(WF.nProtons());
    nn = abs(WF.nNeutrons());
  }

  for ( int i = 0, N = positions.size(); i < N; ++i ) {
    if ( !WF.interNucleonSwing() ) ++isys;

    DipoleStatePtr nucleon;
    tWaveFunctionPtr nwf = WF.nucleon();
    if ( np + nn > 0 ) {
      if ( UseRandom::rndbool(np, nn) ) {
	nwf = WF.proton();
	--np;
      } else {
	nwf = WF.neutron();
	--nn;
      }
    }

    do {
      nucleon = nwf->generate(eh, plus/WF.massNumber());
    } while ( !UseRandom::rndbool(nucleon->weight()) );

    set<PartonPtr> partons;
    for ( vector<DipolePtr>::const_iterator it =
	    nucleon->initialDipoles().begin();
	  it != nucleon->initialDipoles().end();
	  ++it ) {
      partons.insert((**it).partons().first);
      partons.insert((**it).partons().second);
      (**it).colourSystem(isys);
    }

    Parton::Point shift(positions[i].x()/Constants::hbarc,
			positions[i].y()/Constants::hbarc);

    for ( set<PartonPtr>::iterator it = partons.begin();
	  it != partons.end(); ++it )
      (**it).position((**it).position() + shift);

    merge(nucleon);

  }

}

SimpleNucleusState::~SimpleNucleusState() {}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeNoPIOClass<SimpleNucleusState,DIPSY::DipoleState>
  describeDIPSYSimpleNucleusState("DIPSY::SimpleNucleusState",
				  "SimpleNucleus.so");

void SimpleNucleusState::Init() {}

