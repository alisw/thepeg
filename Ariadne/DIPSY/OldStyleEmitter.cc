// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OldStyleEmitter class.
//

#include "OldStyleEmitter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "Parton.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "ThePEG/Utilities/Current.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Repository/UseRandom.h"
#include "gsl/gsl_sf_bessel.h"
#include <limits>

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

OldStyleEmitter::OldStyleEmitter() {}

OldStyleEmitter::~OldStyleEmitter() {}

IBPtr OldStyleEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr OldStyleEmitter::fullclone() const {
  return new_ptr(*this);
}


// Old y generator
 double OldStyleEmitter::oldGenerateY(Dipole dip,double ymin) const {
   double yshift=log(2000.0);
	Energy ppsum = (dip.partons().first->plus() + dip.partons().second->plus())*exp(-yshift);
	return -yshift- log(sqrt(1.0 + sqr(ppsum/(4.0*GeV)))) + 
			sqrt(sqr(log(sqrt(1.0 + sqr(ppsum/(4.0*GeV)))) + ymin+yshift) 
			- log(UseRandom::rnd())/(4.0*alphaBar(Current<DipoleEventHandler>()->rMax())));
 }
 
 //Old xy generator
pair<double,double> OldStyleEmitter::oldGenerateXY(Dipole dip, double y0) const  {
	double ppsum = dip.partons().first->plus()/GeV + dip.partons().second->plus()/GeV;
	double rho = 2.0*exp(-y0)/ppsum ;
  pair<double,double> ret;
  double R1 = UseRandom::rnd();
  double R2 = UseRandom::rnd();
  double ir2 = pow(1.0 + 1.0/(4.0*sqr(rho)),R2) - 1.0;
  double r = 0.0;
  if ( ir2 < std::numeric_limits<double>::epsilon() ) r = rho/sqrt(R2);
  else r = 0.5*sqrt(1.0/ir2); 
  double x = r*cos(2.0*M_PI*R1);
  double y = r*sin(2.0*M_PI*R1);
  ret.first = x;
  if ( UseRandom::rnd() > 0.5 ) ret.first = 1.0 - x;
  ret.second = y;
  return ret;
}



void OldStyleEmitter::generate(Dipole & dip, double ymin, double ymax) const {
	int vetoCounter = 0;

	double y0 = 0.0*log(2000) + max(ymin,max(dip.partons().first->y(),dip.partons().second->y()));
	
   // Generate the r,y and pp values for a possible emission.
 	while ( true ) {
	vetoCounter++;
	//first generate a starting y
	//   cout << "calling generate Y with y0 = " << y0 << endl;
	y0 = oldGenerateY(dip,y0);
	if(thestMode) y0 = 1.0;
	if ( y0 > ymax ) {
	  dip.generatedGluon(new_ptr(Parton()));
		dip.generatedY(y0);
		break;
	}
	
	//OLD GENERATION! SHOULD BE REMOVED EVENTUALLY, BUT KEPT AS BACKUP....
	//generate a xy FIX THIS; WRONG GENERATION!
	pair<double,double> xy = oldGenerateXY(dip,y0);
	
	Energy ppsum = dip.partons().first->plus() + dip.partons().second->plus();
	double pbar = ppsum/4.0/GeV;
	double yshift=log(2000.0);
      y0 += yshift;
      pbar *= exp(-yshift);
    double gfratio = ((log(1.0 + sqr(pbar*exp(y0))))/(log(1.0 + sqr(pbar)) + 2.0*y0))*
      (sqr(xy.first) + sqr(xy.second) + 0.25)*(sqr(1.0 - xy.first) + sqr(xy.second) + 0.25)/
      ((sqr(1.0 - xy.first) + sqr(xy.second))*(sqr(1.0 - xy.first) + sqr(xy.second) + 0.25) +
       (sqr(xy.first) + sqr(xy.second))*(sqr(xy.first) + sqr(xy.second) + 0.25));
         pbar *= exp(yshift);
	 y0 += -yshift;
	if ( gfratio < UseRandom::rnd() )	   
		continue;

	PartonPtr q1 = dip.partons().first;
	PartonPtr q2 = dip.partons().second;
    double rx = xy.first*(q2->position().first*GeV - q1->position().first*GeV) - 
      xy.second*(q2->position().second*GeV - q1->position().second*GeV) + q1->position().first*GeV;
    double ry = xy.first*(q2->position().second*GeV - q1->position().second*GeV) + 
      xy.second*(q2->position().first*GeV - q1->position().first*GeV) + q1->position().second*GeV;  
    double rij = dip.size()*GeV;
    double rin = rij*sqrt(sqr(xy.first) + sqr(xy.second));
    double rjn = rij*sqrt(sqr(1.0 - xy.first) + sqr(xy.second));

    double rmax=Current<DipoleEventHandler>()->rMax()*GeV;
    if(rin+rjn>rij+100) continue;
 if ( sqr(rin*rjn/(rij))*sqr(1.0/rmax)*(sqr(gsl_sf_bessel_K1(rin/rmax)) + 
					sqr(gsl_sf_bessel_K1(rjn/rmax)) - (sqr(rin) + sqr(rjn) - sqr(rij))/(rin*rjn)* gsl_sf_bessel_K1(rin/rmax)*gsl_sf_bessel_K1(rjn/rmax)) < UseRandom::rnd() )
		    continue;

 if ( alphaBar(min(min(rin,rjn),rij)/GeV) < UseRandom::rnd()*alphaBar(rmax/GeV) )
      continue;

    double pt = 2.0/min(rin,rjn);
    double pplus = pt*exp(-y0);
    double pminus = pt*exp(y0);
    if ( q1->plus()/GeV - (rjn/(rjn + rin))*pplus <= 0.0 )
      continue;
    if ( q2->plus()/GeV - (rin/(rin + rjn))*pplus <= 0.0 )
      continue;
    if ( pminus < (2.0/rij)*exp(max(q1->y(),q2->y())) )  //probably wrong, left to reproduce results
      continue;

    double pt1 = max(q1->pT().pt()/GeV, 2.0/rin);
    double pt2 = max(q2->pT().pt()/GeV, 2.0/rjn); 
    if ( log(pt1/(q1->plus()/GeV - (rjn/(rin + rjn))*pplus)) > y0 )
      continue;
    if ( log(pt2/(q2->plus()/GeV - (rin/(rin + rjn))*pplus)) > y0 )
      continue;
	
	
	PartonPtr gluon = new_ptr(Parton());
	gluon->position(Parton::Point(rx/GeV,ry/GeV));
	dip.generatedGluon(gluon);
	dip.generatedY(y0);
	
	break;
	}
}

void OldStyleEmitter::emit(Dipole & d) const {
    TransverseMomentum e1 = TransverseMomentum(1.0*GeV,0.0*GeV);

    double rin = sqrt(d.generatedGluon()->dist2(*d.partons().first))*GeV;
    double rjn = sqrt(d.generatedGluon()->dist2(*d.partons().second))*GeV;
	
    double P1 = rjn/(rin+rjn);
    double P2 = 1.0-P1;
	
    d.generatedGluon()->pT(e1*2.0/min(rin,rjn));
    d.generatedGluon()->y(d.generatedY());
    d.generatedGluon()->plus(d.generatedGluon()->pT().pt()*exp(-d.generatedY()));
    d.generatedGluon()->minus(d.generatedGluon()->pT().pt()*exp(d.generatedY()));
    d.generatedGluon()->oY(d.generatedGluon()->y());
	
    d.partons().first->plus(d.partons().first->plus() - P1*d.generatedGluon()->plus());
    d.partons().first->pT(e1*max(d.partons().first->pT().pt()/GeV, 2.0/rin));
    d.partons().first->y(log(d.partons().first->pT().pt()/d.partons().first->plus()));
    d.partons().first->minus(d.partons().first->pT().pt()*exp(d.partons().first->y()));

    d.partons().second->plus(d.partons().second->plus() - P2*d.generatedGluon()->plus());
    d.partons().second->pT(e1*max(d.partons().second->pT().pt()/GeV, 2.0/rjn));
    d.partons().second->y(log(d.partons().second->pT().pt()/d.partons().second->plus()));
    d.partons().second->minus(d.partons().second->pT().pt()*exp(d.partons().second->y()));

    d.splitDipole(P2);
}



void OldStyleEmitter::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void OldStyleEmitter::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<OldStyleEmitter,DIPSY::Emitter>
  describeDIPSYOldStyleEmitter("DIPSY::OldStyleEmitter", "OldStyleEmitter.so");


void OldStyleEmitter::Init() {

  static ClassDocumentation<OldStyleEmitter> documentation
    ("There is no documentation for the OldStyleEmitter class");

}

