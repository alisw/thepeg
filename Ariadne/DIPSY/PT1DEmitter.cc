// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PT1DEmitter class.
//

#include "PT1DEmitter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "Parton.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "ThePEG/Utilities/Current.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Repository/UseRandom.h"
#include "gsl/gsl_sf_bessel.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

PT1DEmitter::PT1DEmitter() {}

PT1DEmitter::~PT1DEmitter() {}

IBPtr PT1DEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr PT1DEmitter::fullclone() const {
  return new_ptr(*this);
}

void PT1DEmitter::generate(Dipole & dip, double ymin, double ymax) const {
  int vetocounter = 0;
	double y0 = max(ymin,max(dip.partons().first->y(),dip.partons().second->y()));

   // Generate the r,y and pp values for a possible emission.
 	while ( true ) {
	y0 = generateY(& dip,y0,ymax);
	vetocounter++;

	if(thestMode) y0 = 1.0;
	if ( y0 > ymax ) {
	  dip.generatedGluon(new_ptr(Parton()));
		dip.generatedY(y0);
		break;
	}
	Parton::Point p = generateXT(& dip,y0);
	cout << "wrong weight! must take generateXT in account and send it on!!" << endl;
	
	if(OEVeto(& dip,y0,p)) continue;

	 InvEnergy r13 = (dip.partons().first->position() - p).pt();
	 InvEnergy r23 = (dip.partons().second->position() - p).pt();
	double P1 = r23/(r23+r13);
	double P2 = r13/(r13+r23);
	TransverseMomentum v1,v2;

	v1 = 2.0*(p-dip.partons().first->position())/sqr(r13);
	v2 = 2.0*(p-dip.partons().second->position())/sqr(r23);

	Energy pt = max(v1.pt(),v2.pt());
	Energy pplus = pt*exp(-y0);
	Energy pminus = pt*exp(y0);
	if ( dip.partons().first->plus() - P1*pplus <= 0.0*GeV ){
	  // 		cout << "vetoed by limited pp in first parton" << endl;
	  continue;
	}
	if ( dip.partons().second->plus() - P2*pplus <= 0.0*GeV ){
	  // 		cout << "vetoed by limited pp in second parton" << endl;
	  continue;
	}


	bool oldPmOrdering = true;
	if(oldPmOrdering) {
	if ( pminus < 2/dip.size()*exp(max(dip.partons().first->y(),dip.partons().second->y()))) {
	  // 		cout << "vetoed by not p- ordered" << endl;
	  continue;
	}
	}
	else{
	if ( pminus < max(dip.partons().first->minus(),dip.partons().second->minus())) {
	  // 		cout << "vetoed by not p- ordered" << endl;
	  continue;
	}
	}

	Energy pt1 = max(dip.partons().first->pT().pt(),v1.pt());
	Energy pt2 = max(dip.partons().second->pT().pt(),v2.pt());
	Energy plus1 = dip.partons().first->plus() - P1*pplus;
	Energy plus2 = dip.partons().second->plus() - P2*pplus;
	double y1 = log(pt1/plus1);
	double y2 = log(pt2/plus2);
	// *** ATTENTION *** not used:	Energy minus1 = pt1*exp(y1);
	// *** ATTENTION *** not used:  Energy minus2 = pt2*exp(y2);

	//check that also with the recoiled mother dipoles, it is rapidity ordered.
	if ( y1 > y0 ){
	  // 		cout << "vetoed due to recoil -> nonordered y for first parton" << endl;
	  continue;
	}
	if ( y2 > y0 ){
	  // 		cout << "vetoed due to recoil --> nonordered y for second parton" << endl;
	  continue;
	}

// 	cout << "generation passed all vetos, yay!! :) on try number " << vetoCounter << endl;
// 	cout << "r12 = " << r12*GeV << ", Rr13 = " << r13*GeV << ", r23 = " << r23*GeV << ", y0 = " << y0 << endl;
	
	PartonPtr gluon = new_ptr(Parton());
	gluon->position(p);
	dip.generatedGluon(gluon);
	dip.generatedY(y0);
// 	if(vetocounter > 40)
// 	cout << "vetoes: " << vetocounter << endl;
		break;
	}
}

void PT1DEmitter::emit(Dipole & dip) const {
// 	cout << "emitting" << endl;
        Parton::Point p = dip.generatedGluon()->position();
	InvEnergy r13 = (dip.partons().first->position() - p).pt();
	InvEnergy r23 = (dip.partons().second->position() - p).pt();

	TransverseMomentum v1 = 2.0*(p-dip.partons().first->position())/sqr(r13);
	TransverseMomentum v2 = 2.0*(p-dip.partons().second->position())/sqr(r23);
	
	double P1 = r23/(r23+r13);
	double P2 = r13/(r13+r23);
	Energy pt = max(v1.pt(),v2.pt());
	double y = dip.generatedY();

	dip.generatedGluon()->pT(v1/v1.pt()*pt);
	dip.generatedGluon()->y(dip.generatedY());
	dip.generatedGluon()->plus(pt*exp(-y));
	dip.generatedGluon()->minus(pt*exp(y));
	
	dip.partons().first->plus(dip.partons().first->plus() - P1*dip.generatedGluon()->plus());
	if(dip.partons().first->plus() < 0.0*GeV)
	  cout << "omg, negative p1+.. :( " << dip.partons().first->plus()/GeV << endl;
	dip.partons().first->pT(v1/v1.pt()* max(dip.partons().first->pT().pt(),v1.pt()));
	dip.partons().first->y(log(dip.partons().first->pT().pt()/dip.partons().first->plus()));
	dip.partons().first->minus(dip.partons().first->pT().pt()*exp(dip.partons().first->y()));
   
	dip.partons().second->pT(v2/v2.pt()*max(dip.partons().second->pT().pt(),v2.pt()));
	dip.partons().second->plus(dip.partons().second->plus() - P2*dip.generatedGluon()->plus());
	if(dip.partons().second->plus() < 0.0*GeV) cout << "omg, negative p2+.. :(" << endl;
	dip.partons().second->y(log(dip.partons().second->pT().pt()/dip.partons().second->plus()));
	dip.partons().second->minus(dip.partons().second->pT().pt()*exp(dip.partons().second->y()));
	
	// 	cout << "new p+ sum is " << (d.partons().first->plus()+d.partons().second->plus()+d.generatedGluon()->plus())/GeV << endl;
	// 	cout << "y of the new gluon is " << d.generatedGluon()->y() << endl;
	
	// *** ATTENTION *** Implement this.
}



void PT1DEmitter::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void PT1DEmitter::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<PT1DEmitter,DIPSY::Emitter>
  describeDIPSYPT1DEmitter("DIPSY::PT1DEmitter", "PT1DEmitter.so");

void PT1DEmitter::Init() {

  static ClassDocumentation<PT1DEmitter> documentation
    ("There is no documentation for the PT1DEmitter class");

}

