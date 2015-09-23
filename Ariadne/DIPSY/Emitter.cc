// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Emitter class.
//

//Commented

#include "Emitter.h"
#include "Parton.h"
#include "ShadowParton.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "EffectiveParton.h"
#include "ThePEG/Utilities/Current.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DebugItem.h"

#include <iostream>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Emitter.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "gsl/gsl_sf_bessel.h"

using namespace DIPSY;

Emitter::~Emitter() {}

/*
 * Just a shortcut to access the input parameter from the Current handler.
 */
InvEnergy Emitter::rMax() const {
  return theRMax > 0.0*InvGeV? theRMax: Current<DipoleEventHandler>()->rMax();
}

/*
 * Shortcut to access running alpha from the Current handler.
 */
double Emitter::alphaS(InvEnergy r) const {
  return Current<DipoleEventHandler>()->alphaS(r);
}

/*
 * The veto check for the overestimates done previously.
 * This calculates the overestimated emission probability, and checks the
 * actual probability it should have been in this point, and then returns
 * true or false randomly depending on the ratio.
 * Note that this is not the final distribution, as there are still more
 * vetos done later.
 */
bool Emitter::OEVeto(DipolePtr dip, double y0, Parton::Point p) const {
  //set up dipole sizes.
  InvEnergy R13 = (dip->effectivePartons().first->position() - p).pt();
  InvEnergy R23 = (dip->effectivePartons().second->position() - p).pt();
  InvEnergy R12 = dip->size();

  //Too large dipoles messes up the bessel functions, and will not make it
  //through the p- veto anyways, so can be vetoed already now, to avoid
  //errors from the bessel function.
  if(R13*GeV > 100.0) return true;
  if(R23*GeV > 100.0) return true;
  if(R12*GeV > 100.0) return false;

  //some more set up.
  double bess13 = gsl_sf_bessel_K1(R13/rMax());
  double bess23 = gsl_sf_bessel_K1(R23/rMax());
  
  //this is the correct probability distribution used in the paper.
  Energy2 correctDist = alphaBar(min(min(R13,R23),R12))/(M_PI*2.0*sqr(rMax()))*
    (sqr(bess13) + sqr(bess23) - (sqr(R13) + sqr(R23) - sqr(R12))/
     (R13*R23)*bess13*bess23);

  // If we create dipoles which are larger than the original, the
  // original partons will be distributed as dpt2/pt4. Here we may
  // include an extra fudge factor to emulate a matrix element
  // correction.
  if ( R13 > R12 && R23 > R12 && Current<DipoleEventHandler>()->fudgeME() )
    correctDist *= 1.0 - 1.0/(1.0 + cosh(dip->partons().first->y() -
					 dip->partons().second->y()));

  //Now calculate the overestimated distribution in the same way as it is done in gerateY() and generateXT().
  double Coe = alphaBar(R12/2.0)/M_PI*
    sqr(gsl_sf_bessel_K1(R12/(2.0*rMax())))*
    sqr(R12/(2.0*rMax()))*
    (R12/(2.0*theRScale)+2.0);
  if(R12>2.0*rMax()) //double check that this is ok!!!!
    Coe *= sqrt(R12/(2.0*rMax()))*
      exp(R12/(rMax())-2.0);

  //this is the overestime
  Energy2 overestimate = Coe*(1.0/(sqr(R13)*(R13/theRScale+2.0)) + 
                                   1.0/(sqr(R23)*(R23/theRScale+2.0)));

  //if the overestimate is ok, this should never happen.
  if( correctDist/overestimate > 1.0 )
    Throw<EmitterException>()
      << "In DIPSY::Emitter " << name() << ": XT keep prob is " << correctDist/overestimate 
      << " (r12 = " << R12*GeV << ", Rr13 = " << R13*GeV << ", r23 = " << R23*GeV << ")"
      << Exception::warning;

  //generate a random number and check if veto or not.
  if(correctDist/overestimate>UseRandom::rnd())
    return false;
  else
    return true;
}

/*
 * The veto check for the overestimates done previously using shadows.
 * This calculates the overestimated emission probability, and checks the
 * actual probability it should have been in this point, and then returns
 * true or false randomly depending on the ratio.
 * Note that this is not the final distribution, as there are still more
 * vetos done later.
 */
bool Emitter::OEVeto(DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2,
		     double y0, Parton::Point p) const {
  //set up dipole sizes.
  InvEnergy R13 = (dip->partons().first->position() - p).pt();
  InvEnergy R23 = (dip->partons().second->position() - p).pt();
  InvEnergy R12 = dip->size();

  //Too large dipoles messes up the bessel functions, and will not make it
  //through the p- veto anyways, so can be vetoed already now, to avoid
  //errors from the bessel function.
  if ( R13*GeV > 100.0 || R23*GeV > 100.0 ) return true;
  if ( R12*GeV > 100.0 ) return false;

  //some more set up.
  double bess13 = gsl_sf_bessel_K1(R13/rMax());
  double bess23 = gsl_sf_bessel_K1(R23/rMax());
  
  //this is the correct probability distribution used in the paper.
  Energy2 correctDist = alphaBar(min(min(R13,R23),R12))/(M_PI*2.0*sqr(rMax()))*
    (sqr(bess13) + sqr(bess23) - (sqr(R13) + sqr(R23) - sqr(R12))/
     (R13*R23)*bess13*bess23);

  // If we create dipoles which are larger than the original, the
  // original partons will be distributed as dpt2/pt4. Here we may
  // include an extra fudge factor to emulate a matrix element
  // correction.
  if ( R13 > R12 && R23 > R12 && Current<DipoleEventHandler>()->fudgeME() )
    correctDist *= 1.0 - 1.0/(1.0 + cosh(dip->partons().first->y() -
					 dip->partons().second->y()));

  //Now calculate the overestimated distribution in the same way as it is done in gerateY() and generateXT().
  double Coe = alphaBar(R12/2.0)/M_PI*
    sqr(gsl_sf_bessel_K1(R12/(2.0*rMax())))*
    sqr(R12/(2.0*rMax()))*
    (R12/(2.0*theRScale)+2.0);
  if(R12>2.0*rMax()) //double check that this is ok!!!!
    Coe *= sqrt(R12/(2.0*rMax()))*
      exp(R12/(rMax())-2.0);

  //this is the overestime
  Energy2 overestimate = Coe*(1.0/(sqr(R13)*(R13/theRScale+2.0)) + 
                                   1.0/(sqr(R23)*(R23/theRScale+2.0)));

  //if the overestimate is ok, this should never happen.
  //TODO: write a proper error message.
  if ( correctDist/overestimate > 1.0 )
    Throw<EmitterException>()
      << "In DIPSY::Emitter " << name() << ": XT keep prob is " << correctDist/overestimate 
      << " (r12 = " << R12*GeV << ", Rr13 = " << R13*GeV << ", r23 = " << R23*GeV << ")"
      << Exception::warning;

  //generate a random number and check if veto or not.
  return !UseRandom::rndbool(correctDist/overestimate);

}

/*
 * This corrects the overestimated emission probability in Y used in generateY().
 * The correct emission rate is calculated for the rapidity in question and
 * compared to the overestimate rate used in generateY().
 * Note that this is not the final correct y-distribution, as there are still several
 * vetos to be checked, it is just a step closer to the correct distribution.
 */
bool Emitter::YVeto(double y,DipolePtr dip, double Coe, double rateOE) const {
  tEffectivePartonPtr p1 = dip->effectivePartons().first;
  tEffectivePartonPtr p2 = dip->effectivePartons().second;

  //calculate the correct rate at the rapidity y in question.
  //TODO: member calculating cutoff, interface before/after recoil.
  //Before recoil corresponds to order wrt propagator, after wrt final gluon.
  InvEnergy cutoff = 0.99*0.5*
    min(pTScale()/(p1->plus()*exp(y) + p1->pT().pt()/2.0),
    min(pTScale()/(p2->plus()*exp(y) + p2->pT().pt()/2.0),
        dip->size()/4.0)); //if ordered AFTER recoil
  // InvEnergy cutoff = 0.99*
  //   min(2.0/3.0*pTScale()*exp(-y)/p1->plus(),
  //   min(2.0/3.0*pTScale()*exp(-y)/p2->plus(),
  //       dip->size()/4.0)); //if ordered BEFORE recoil

  //calculate overestimated rate according the what is used in generateY().
  double rate = 2.0*M_PI*Coe*log(1.0 + 2.0*theRScale/cutoff);

  //generate rnd() and veto.
  if(rate/rateOE > UseRandom::rnd())		//if exp(-...) big enough
    return false; 				//keep
  else
    return true;					//else veto
}

/*
 * This corrects the overestimated emission probability in Y used in generateY().
 * The correct emission rate is calculated for the rapidity in question and
 * compared to the overestimate rate used in generateY().
 * Note that this is not the final correct y-distribution, as there are still several
 * vetos to be checked, it is just a step closer to the correct distribution.
 */
bool Emitter::YVeto(double y, DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2,
		    double Coe, double rateOE) const {
  //calculate the correct rate at the rapidity y in question.
  //TODO: member calculating cutoff, interface before/after recoil.
  //Before recoil corresponds to order wrt propagator, after wrt final gluon.
  InvEnergy cutoff = 0.99*0.5*
    min(pTScale()/(sp1->plus0()*exp(y) + sp1->pT0().pt()/2.0),
    min(pTScale()/(sp2->plus0()*exp(y) + sp2->pT0().pt()/2.0),
        dip->size()/4.0));

  //calculate overestimated rate according the what is used in generateY().
  double rate = 2.0*M_PI*Coe*log(1.0 + 2.0*theRScale/cutoff);

  return !UseRandom::rndbool(rate/rateOE);

}

/*
 * Generates a rapidity for the next emission.
 * The distribution f(y) in this member is using is an overestimate, and will be
 * vetoed down to correct distribution in the checks of the x_T distribution.
 *
 * Already the overestimate f(y) is too messy to generate directly though,
 * so it has to be done through a further overestimate g(y) > f(y) which is
 * then vetoed (only depending on y) down to f(y) with YVeto().
 * Look at the writeup for more details on f() and g().
 */
double Emitter::
generateY(DipolePtr dip, double ymin, double ymax) const {
  //set up shorter names.
  tEffectivePartonPtr p1 = dip->effectivePartons().first;
  tEffectivePartonPtr p2 = dip->effectivePartons().second;
  InvEnergy r = dip->size();

  //To large dipoles will mess up the bessel function, so just shut them
  // off by emitting after the maximum y.
  if ( r*GeV > 100.0 ) return ymax + 1.0;

  //Calculate the cutoff in transverse distance. See writeup for details.
  InvEnergy cutoff = 0.99*0.5*
    min(pTScale()/(p1->plus()*exp(ymin) + p1->pT().pt()/2.0),
    min(pTScale()/(p2->plus()*exp(ymin) + p2->pT().pt()/2.0),
        dip->size()/4.0)); //if ordered AFTER recoil
  // InvEnergy cutoff = 0.99*
  //   min(2.0/3.0*pTScale()*exp(-ymin)/p1->plus(),
  //   min(2.0/3.0*pTScale()*exp(-ymin)/p2->plus(),
  //       dip->size()/4.0)); //if ordered BEFORE recoil

  //Calculate the overestimated coefficient. See writeup for details.
  double Coe = alphaBar(r/2.0)/M_PI*
    sqr(gsl_sf_bessel_K1(r/(2.0*rMax())))*
    sqr(r/(2.0*rMax()))*
    (r/(2.0*theRScale)+2.0);

  //For assymptotically large dipoles, the overestimated emission rate grows
  //faster than Coe, so add an extra factor for large dipoles.
  //this will have to be removed in the veto.
  if(r>2.0*rMax()) //double check that this is ok!!!
    Coe *= sqrt(r/(2.0*rMax()))*
      exp(r/(rMax())-2.0);

  //calculate the overestimated rate of emission (per unit rapidity).
  double rateOE = 2.0*M_PI*Coe*log(1.0 + 2.0*theRScale/cutoff);

  //generate the rapidity. (see writeup for details.)
  double y = ymin + log(1.0-log(UseRandom::rnd())/rateOE);

  //veto to get to correct distribution, recur if vetoed.
  if(YVeto(y,dip,Coe,rateOE*exp(y-ymin)) && y < ymax)
    return generateY(dip, y, ymax);
  else{
    return y;
  }
 }

/*
 * Generates a rapidity for the next emission using shadows.
 * The distribution f(y) in this member is using is an overestimate, and will be
 * vetoed down to correct distribution in the checks of the x_T distribution.
 *
 * Already the overestimate f(y) is too messy to generate directly though,
 * so it has to be done through a further overestimate g(y) > f(y) which is
 * then vetoed (only depending on y) down to f(y) with YVeto().
 * Look at the writeup for more details on f() and g().
 */
double Emitter::generateY(DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2,
			  double ymin, double ymax) const {
  //set up shorter names.
  InvEnergy r = dip->size();

  //To large dipoles will mess up the bessel function, so just shut them
  // off by emitting after the maximum y.
  if ( r*GeV > 100.0 ) return ymax + 1.0;

  //Calculate the cutoff in transverse distance. See writeup for details.
  InvEnergy cutoff = 0.99*0.5*
    min(pTScale()/(sp1->plus0()*exp(ymin) + sp1->pT0().pt()/2.0),
	min(pTScale()/(sp2->plus0()*exp(ymin) + sp2->pT0().pt()/2.0), r/4.0));

  //Calculate the overestimated coefficient. See writeup for details.
  double Coe = alphaBar(r/2.0)/M_PI*
    sqr(gsl_sf_bessel_K1(r/(2.0*rMax())))*
    sqr(r/(2.0*rMax()))*
    (r/(2.0*theRScale)+2.0);

  //For assymptotically large dipoles, the overestimated emission rate grows
  //faster than Coe, so add an extra factor for large dipoles.
  //this will have to be removed in the veto.
  if(r>2.0*rMax()) //double check that this is ok!!!
    Coe *= sqrt(r/(2.0*rMax()))*
      exp(r/(rMax())-2.0);

  //calculate the overestimated rate of emission (per unit rapidity).
  double rateOE = 2.0*M_PI*Coe*log(1.0 + 2.0*theRScale/cutoff);

  //generate the rapidity. (see writeup for details.)
  double y = ymin + log(1.0-log(UseRandom::rnd())/rateOE);

  //veto to get to correct distribution, recur if vetoed.
  if( YVeto(y, dip, sp1, sp2, Coe, rateOE*exp(y-ymin)) && y < ymax )
    return generateY(dip, sp1, sp2, y, ymax);

  return y;

 }

/*
 * Generate the transverse position. First generate the distance r
 * larger than the cutoff used in generateY(), then decide around which
 * parton, and then the angle. Distributions found in writeup.
 */
Parton::Point Emitter::generateXT(DipolePtr dip, double y0) const {
  //set up
  tEffectivePartonPtr p1 = dip->effectivePartons().first;
  tEffectivePartonPtr p2 = dip->effectivePartons().second;

  //Calculate the cutoff.
  //TODO: make cutoff member, taking dipole and y as arguments.
  InvEnergy cutoff = 0.99*0.5*
    min(pTScale()/(p1->plus()*exp(y0) + p1->pT().pt()/2.0),
    min(pTScale()/(p2->plus()*exp(y0) + p2->pT().pt()/2.0),
        dip->size()/4.0)); //if ordered AFTER recoil
  // InvEnergy cutoff = 0.99*
  //   min(2.0/3.0*pTScale()*exp(-y0)/p1->plus(),
  //   min(2.0/3.0*pTScale()*exp(-y0)/p2->plus(),
  //       dip->size()/4.0)); //if ordered BEFORE recoil

  //The normalisation. See writeup for details.
  double C2 = 2.0/log(1.0+2.0*theRScale/cutoff);	//normalises g(z)

  //Generate the distance.
  InvEnergy r = 2.0*theRScale/(exp(2.0*UseRandom::rnd()/C2)-1.0);
  //tested to give correct distribution C2/(r^3/theRScale+2r^2)

  //generate the angle.
  double phi = 2.0*M_PI*UseRandom::rnd();

  //decide which parton, and create the Point to be returned.
  if(UseRandom::rnd()>0.5)
    return Parton::Point(p1->position().first + cos(phi)*r,
				   p1->position().second + sin(phi)*r);
  else
    return Parton::Point(p2->position().first + cos(phi)*r,
				   p2->position().second + sin(phi)*r);
}

/*
 * Generate the transverse position. First generate the distance r
 * larger than the cutoff used in generateY(), then decide around which
 * parton, and then the angle. Distributions found in writeup.
 */
Parton::Point
Emitter::generateXT(DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2, double y0) const {

  //Calculate the cutoff.
  //TODO: make cutoff member, taking dipole and y as arguments.
  InvEnergy cutoff = 0.99*0.5*
    min(pTScale()/(sp1->plus0()*exp(y0) + sp1->pT0().pt()/2.0),
    min(pTScale()/(sp2->plus0()*exp(y0) + sp2->pT0().pt()/2.0),
        dip->size()/4.0)); //if ordered AFTER recoil

  //The normalisation. See writeup for details.
  double C2 = 2.0/log(1.0 + 2.0*theRScale/cutoff);	//normalises g(z)

  //Generate the distance.
  InvEnergy r = 2.0*theRScale/(exp(2.0*UseRandom::rnd()/C2)-1.0);
  //tested to give correct distribution C2/(r^3/theRScale+2r^2)

  //generate the angle.
  double phi = 2.0*M_PI*UseRandom::rnd();

  //decide which parton, and create the Point to be returned.
  tPartonPtr pp =
    UseRandom::rndbool()? dip->partons().first: dip->partons().second;
    return Point(pp->position().first + cos(phi)*r,
		 pp->position().second + sin(phi)*r);

}

/*
* the main function to be called from outside. Other functions are mainly to help
* this one. To get the emission right, several layers of overestimates and vetoes
* are used.
*/
void Emitter::generate(Dipole & dip, double ymin, double ymax) const {

  if ( dip.partons().first->shadow() ) {
    generateWithShadows(dip, ymin, ymax);
    return;
  }

  //Lowest allowed emission rapidity. Not before any of the partons
  //in the dipole, and not before the supplied ymin.
  double y0 = max(ymin,max(dip.partons().first->y(),dip.partons().second->y()));

  //define the resolved effective partons, as function of the dipole size.
  //Later, when the transverse position is decided, this will have to
  //be recalculated with the actual resolution ranges.
  //This is an overestimate, as later lowering the range can only
  //decrease the allowed phase space through less p+, and increased p_T.
  dip.effectivePartons(EffectiveParton::create(*(dip.partons().first),
					       dip.size()/2.0),
		       EffectiveParton::create(*(dip.partons().second),
					       dip.size()/2.0));

  //Some renaming to save space later.
  tEffectivePartonPtr p1 = dip.effectivePartons().first;
  tEffectivePartonPtr p2 = dip.effectivePartons().second;

  //We will go through the while loop until something passes all the vetoes
  //or pass the maximum allowed rapidity ymax. However, put a limit
  //on the number of trials not to get stuck. So define the count to keep
  //track. Possibly a for loop would be more appropriate at this point...
  while ( true ) {
    //reset the range to the overestimate, as we do not know
    //what x_T will eb chosen this trial.
    p1->setRange( dip.size()/2.0 );
    p2->setRange( dip.size()/2.0 );

    //generate a rapidity, according to an overestimate.
    y0 = generateY(& dip, y0, ymax);


    //if the generated rapidity is above the limit, return an empty pointer
    if ( y0 > ymax ) {
      dip.generatedGluon(new_ptr(Parton()));
      dip.generatedY(y0);
      break;
    }

    //generate a transverse position, using an overestimated
    //phase space.
    Parton::Point p = generateXT(& dip, y0);

    //Bring down the overestimate in the distribution
    //to get the right amplitude. After this, only phase space limitations
    //left.
    if(OEVeto(& dip,y0,p)) {
      continue;
    }

    //Set up notation.
    InvEnergy r13 = (p1->position() - p).pt();
    InvEnergy r23 = (p2->position() - p).pt();
    InvEnergy r12 = dip.size();

    //this is to what extent the first parent emitted, and to what extent the
    //second one. Here, it is some kind of superposition where they are emitting
    //together, so both parents recoil a bit, etc.
    double P1 = sqr(r23)/(sqr(r23) + sqr(r13));
    double P2 = sqr(r13)/(sqr(r13) + sqr(r23));

    //check how much is actually resolved, small emissions resolve more
    if ( Current<DipoleEventHandler>()->emitter().rangeMode() != 1 ) {
      if ( r13 < r12/2.0 )
	p1->setRange( r13 );
      if ( r23 < r12/2.0 )
	p2->setRange( r23 );
    }

    //The transverse momentum of the new parton from the parents.
    TransverseMomentum v1 = pTScale()*(p - p1->position())/sqr(r13);
    TransverseMomentum v2 = pTScale()*(p - p2->position())/sqr(r23);

    //calculate the 4-momentum of the emitted gluon.
    TransverseMomentum pt = v1 + v2;
    Energy pplus = pt.pt()*exp(-y0);
    Energy pminus = pt.pt2()/pplus;

    //check that there's enough p+
    if ( p1->plus() - P1*pplus <= 0.0*GeV ){
      continue;
    }
    if ( p2->plus() - P2*pplus <= 0.0*GeV ){
      continue;
    }

    //calculate the new 4-momenta of the recoiled parents.
    TransverseMomentum pt1 = p1->pT() - v1;
    TransverseMomentum pt2 = p2->pT() - v2;
    Energy plus1 = p1->plus() - P1*pplus;
    Energy plus2 = p2->plus() - P2*pplus;

    double y1 = log(pt1.pt()/plus1);
    double y2 = log(pt2.pt()/plus2);
    Energy minus1 = pt1.pt2()/plus1;
    Energy minus2 = pt2.pt2()/plus2;

    //this option use the p- ordering of the non-effective parton,
    //so p- should be calculated from the non-effective parton.
    if ( theMinusOrderingMode == 1 ) {
      minus1 = pt1.pt()*exp(dip.partons().first->y());
      minus2 = pt2.pt()*exp(dip.partons().second->y());
    }

    //in the first option, p- is required to be fully ordered with both parents.
    //in the second option, the parton has to be ordered mainly to the parton that
    //emitted it most, so it is allowed to be unordered by a factor P1 and P2.
    //PSInflation is a plain factor that the emission is allowed to be unordered by.
    if ( bothOrderedEvo() ) {
      if ( pminus*thePSInflation < max(minus1, minus2)*thePMinusOrdering ) continue;
    }
    else {
      if ( pminus*thePSInflation < max(P1*minus1, P2*minus2)*thePMinusOrdering ) continue;
    }
    
    //check rapidity ordered with the two recoiled parents.
    if ( y1 > y0 ) {
      continue;
    }
    if ( y2 > y0 ) {
      continue;
    }

    //check that none of the internal partons recoiled over ymax.
    if ( p1->recoilsOverYMax( v1, P1*pplus, ymax ) ) {
      continue;
    }
    if ( p2->recoilsOverYMax( v2, P2*pplus, ymax ) ) {
      continue;
    }

    //This is a self-consistency check.
    //If things are donw correctly, then emissions at the limit of the
    //cutoff in x_T used in generateXT() and generateY() should all be vetoed.
    //if emissions at the limit go through, it means that we are missing part of
    //the phase space just within the cutoff, which is bad.
    //
    //to check this, the cutoff is decreased by a factor 0.99 in generateY() etc above,
    //to create the occasional emission in a region that should always be vetoed if the
    //overestimates are fair. This is a check at the end so that none of these emission
    //between 1 and 0.99 of the cutoff went through.
    //TODO: match this with cutoff above
    //TODO: proper error message.
    if ( min(r13, r23) < min(2.0/3.0*pTScale()/(p1->plus()*exp(y0)),
			 min(2.0/3.0*pTScale()/(p2->plus()*exp(y0)),
			     dip.size()/4.0)) )
      Throw<EmitterException>()
	<< "In DIPSY::Emitter " << name() << ": Emission below r-cutoff passed vetos!"
	<< " (r12 = " << ounit(r12, InvGeV) << ", r13 = " << ounit(r13, InvGeV)
	<< ", r23 = " << ounit(r23, InvGeV)
	<< "v1/pt = " << double(min(v1.pt(), v2.pt())/pt.pt()) << ")"
	<< Exception::warning;
    
    // passed. Set up a new parton for the dipole with the information of the emission.
    PartonPtr gluon = new_ptr(Parton());
    gluon->position(p);
    dip.generatedGluon(gluon);
    dip.generatedY(y0);

    //leave the while loop.
    break;
  }
}

/*
* the main function to be called from outside when using
* shadows. Other functions are mainly to help this one. To get the
* emission right, several layers of overestimates and vetoes are used.
*/
void Emitter::generateWithShadows(Dipole & dip, double ymin, double ymax) const {

  static DebugItem trace("DIPSY::Trace", 9);

  if ( trace ) cerr << "Gen  " << dip.tag() << endl;

  // Handy pointers to the partons.
  tPartonPtr p1 = dip.partons().first;
  tPartonPtr p2 = dip.partons().second;

  //First check how far down the history we must consider previous shadows.
  tSPartonPtr sp1 =
    dip.partons().first->shadow()->resolve(dip.size2()/4.0, dip.partons().second);
  tSPartonPtr sp2 =
    dip.partons().second->shadow()->resolve(dip.size2()/4.0, dip.partons().first);

  //Lowest allowed emission rapidity. Not before any of the partons
  //in the dipole, and not before the supplied ymin.
  double y0 = max(ymin, max(p1->y(), p2->y()));

  //We will go through the while loop until something passes all the vetoes
  //or pass the maximum allowed rapidity ymax. However, put a limit
  //on the number of trials not to get stuck. So define the count to keep
  //track. Possibly a for loop would be more appropriate at this point...
  while ( true ) {
    //generate a rapidity, according to an overestimate.
    y0 = generateY(&dip, sp1, sp2, y0, ymax);


    //if the generated rapidity is above the limit, return an empty pointer
    if ( y0 > ymax ) {
      dip.generatedGluon(new_ptr(Parton()));
      dip.generatedY(y0);
      return;
    }

    //generate a transverse position, using an overestimated
    //phase space.
    Parton::Point p = generateXT(&dip, sp1, sp2, y0);

    //Bring down the overestimate in the distribution
    //to get the right amplitude. After this, only phase space limitations
    //left.
    if ( OEVeto(&dip, sp1, sp2, y0, p) ) {
      continue;
    }

    //Set up notation.
    InvEnergy2 r13 = (p1->position() - p).pt2();
    InvEnergy2 r23 = (p2->position() - p).pt2();
    InvEnergy2 r12 = dip.size2();

    // We need to already here decide which parton is emitting.
    bool firstsplit = UseRandom::rndbool(r23/(r23 + r13));
    tPartonPtr pe = firstsplit? p1: p2;
    tPartonPtr pr = firstsplit? p2: p1;
    tSPartonPtr spe = firstsplit? sp1: sp2;
    InvEnergy2 re3 = firstsplit? r13: r23;

    // The emitting parton may be further resolved by the emission.
    tSPartonPtr spre =  re3 < r12/4.0? pe->shadow()->resolve(re3, pr): spe;

    //The transverse momentum of the new parton from the emitter.
    TransverseMomentum pt = pTScale()*(p - pe->position())/re3;
    Energy pplus = pt.pt()*exp(-y0);
    Energy pminus = pt.pt2()/pplus;

    //check that there's enough p+
    if ( pplus > spre->plus0() ) continue;

    //calculate the new 4-momenta of the recoiled parents.
    TransverseMomentum pte = spre->pT() - pt;
    Energy pluse = spre->plus0() - pplus;

    double ye = log(pte.pt()/pluse);
    Energy minuse = pte.pt2()/pluse;
    if ( theMinusOrderingMode >= 2 ) {
      ShadowParton::Propagator prop =
	spre->propagator(min(re3, r12/4.0), pr, -1);
      if ( prop.fail ) continue;
      pluse = prop.p.plus() - pplus;
      if ( pluse <= ZERO ) continue;
      pte = TransverseMomentum(prop.p) - pt;
      ye = log(pte.pt()/pluse);
      minuse = pte.pt2()/pluse;
      if ( theMinusOrderingMode >= 3 ) {
	if ( firstsplit ) {
	  if ( pplus > prop.colpos || pminus < prop.colneg ) continue;
	} else {
	  if ( pplus > prop.acopos || pminus < prop.aconeg ) continue;
	}
      }

      // *** TODO *** fix masses!
    }
    else if ( theMinusOrderingMode == 1 )
      minuse = pte.pt()*exp(pe->y());

    // Ensure approximate minus ordering.
    if ( pminus*thePSInflation < minuse*thePMinusOrdering ) continue;

    //check rapidity ordered with the two recoiled parents.
    if ( ye > y0 ) continue;

    // passed. Set up a new parton for the dipole with the information of the emission.
    PartonPtr gluon = new_ptr(Parton());
    gluon->position(p);
    gluon->mainParent(pe);
    dip.generatedGluon(gluon);
    dip.generatedY(y0);

    // And we're done!
    return;
  }
}

/*
 * Maked the dipole actually emit the dipole prepared in generate(), setting up
 * parents, children, changing 4-momenta etc.
 */
void Emitter::emit(Dipole & d) const {
  //some notation.
  if ( d.partons().first->shadow() ) {
    emitWithShadows(d);
    return;
  }
  PartonPtr p = d.generatedGluon();
  InvEnergy r13 = (d.partons().first->position() - p->position()).pt();
  InvEnergy r23 = (d.partons().second->position() - p->position()).pt();

  //The recoils.
  //TODO: use recoil().
  TransverseMomentum v1 = pTScale()*(p->position() - d.partons().first->position())/sqr(r13);
  TransverseMomentum v2 = pTScale()*(p->position() - d.partons().second->position())/sqr(r23);
  
  //the ratio each parent emitted the gluon with.
  double P1 = sqr(r23)/(sqr(r13) + sqr(r23));
  double P2 = 1.0 - P1;

  //set 4-momenta for emission.
  p->pT(v1 + v2);
  p->y(d.generatedY());
  p->plus(p->pT().pt()*exp(-d.generatedY()));
  p->minus(p->pT().pt()*exp(d.generatedY()));
  p->oY(p->y());

  //change the dipole structure, generate new colour etc.
  d.splitDipole(P2);

  //perform the recoil on the effective partons.
  //the resolution ranges should be set since the generation of the emission.
  d.effectivePartons().first->recoil( v1, P1*p->plus() );
  d.effectivePartons().second->recoil( v2, P2*p->plus() );
}

void Emitter::emitWithShadows(Dipole & d) const {
  static DebugItem trace("DIPSY::Trace", 9);

  if ( trace ) cerr << "Emit " << d.tag();

  //some notation.
  PartonPtr p = d.generatedGluon();

  // Decide which parton is the emitter.
  bool em1 = ( p->mainParent() == d.partons().first );
  tPartonPtr emitter = em1? d.partons().first: d.partons().second;
  tPartonPtr recoiler = em1? d.partons().second: d.partons().first;
  InvEnergy2 res = emitter->dist2(*p);
  tSPartonPtr sem = emitter->shadow()->resolve(min(res, d.size2()/4.0), recoiler);

  TransverseMomentum rec =
    pTScale()*(p->position() - emitter->position())/res;

  //change the dipole structure, generate new colour etc.
  d.splitDipole(em1? 0.0: 1.0);

  //set 4-momenta for emission.
  p->pT(rec);
  p->y(d.generatedY());
  p->plus(p->pT().pt()*exp(-d.generatedY()));
  p->minus(p->pT().pt()*exp(d.generatedY()));
  p->oY(p->y());

  if ( theMinusOrderingMode >= 2 ) {
    ShadowParton::Propagator prop =
      sem->propagator(min(res, d.size2()/4.0), recoiler, -1);
    emitter->pT(TransverseMomentum(prop.p) - rec);
    emitter->plus(prop.p.plus()  - p->plus());
    emitter->y(log(emitter->mt()/emitter->plus()));
    emitter->minus(emitter->mt2()/emitter->plus());
  } else {
    emitter->pT(sem->pT() - rec);
    emitter->plus(sem->plus() - p->plus());
    emitter->y(log(emitter->pT().pt()/emitter->plus()));
    emitter->minus(emitter->mt2()/emitter->plus());
  }

  emitter->shadow()->setupEmission(*emitter, *p, *recoiler);

  if ( trace ) cerr << " -> " << d.children().first->tag()
		    << " + " << d.children().second->tag() << endl;

}

void Emitter::persistentOutput(PersistentOStream & os) const {
  os << ounit(thePSInflation, 1.0)
     << ounit(thePMinusOrdering, 1.0)
     << ounit(theRScale, InvGeV) 
     << ounit(theRMax, InvGeV) 
     << ounit(thePTScale, 1.0)
     << ounit(theBothOrderedEvo, true)
     << ounit(theBothOrderedInt, true)
     << ounit(theBothOrderedFS, true)
     << ounit(theRangeMode, 1)
     << ounit(theMinusOrderingMode, 1);
}

void Emitter::persistentInput(PersistentIStream & is, int) {
  is >> iunit(thePSInflation, 1.0)
     >> iunit(thePMinusOrdering, 1.0)
     >> iunit(theRScale, InvGeV) 
     >> iunit(theRMax, InvGeV)
     >> iunit(thePTScale, 1.0)
     >> iunit(theBothOrderedEvo, true)
     >> iunit(theBothOrderedInt, true)
     >> iunit(theBothOrderedFS, true)
     >> iunit(theRangeMode, 1)
     >> iunit(theMinusOrderingMode, 1);
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<Emitter,HandlerBase> 
describeDIPSYEmitter("DIPSY::Emitter", "libAriadne5.so libDIPSY.so");

void Emitter::Init() {

  static ClassDocumentation<Emitter> documentation
    ("The Emitter class is responsible for generating and performing "
     "emissions from dipoles. This base class does the default emission "
     "strategy.");


  static Parameter<Emitter,InvEnergy> interfaceRScale
    ("RScale",
     "Constant used in overestimate of emission probability.",
     &Emitter::theRScale, InvGeV, 1.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<Emitter,double> interfacePSInflation
    ("PSInflation",
     "How much p+- ordering should be enforced. 1 is normal ordering, high numbers "
     "are no ordering(energy must still be conserved though), "
     "low number demands stronger ordering",
     &Emitter::thePSInflation, 1.0, 1.0, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Parameter<Emitter,double> interfacePMinusOrdering
    ("PMinusOrdering",
     "An extra factor strengthening the p- ordering on top of the PSInflation."
     "A large value will suppress large dipoles more.",
     &Emitter::thePMinusOrdering, 1.0, 1.0, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Parameter<Emitter,double> interfacePTScale
    ("PTScale",
     "If pT is 1/r or 2/r.",
     &Emitter::thePTScale, 1.0, 1.0, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Parameter<Emitter,InvEnergy> interfaceRMax
    ("RMax",
     "The confinement scale (in iverse GeV). If set to zero, "
     "the value of <interface>DipoleEventHandler::RMax</interface> of the "
     "controlling event handler will be used.",
     &Emitter::theRMax, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0*InvGeV,
     true, false, Interface::lowerlim);

  static Switch<Emitter,bool> interfaceBothOrderedEvo
    ("BothOrderedEvo",
     "If an emission should be fully ordered with both its parents."
     "otherwise it will be weighted according to how close they are.",
     &Emitter::theBothOrderedEvo, false, true, false);
  static SwitchOption interfaceBothOrderedEvoTrue
    (interfaceBothOrderedEvo,"True","both parents fully ordered.",true);
  static SwitchOption interfaceBothOrderedEvoFalse
    (interfaceBothOrderedEvo,"False","weighted by distance.",false);

  static Switch<Emitter,bool> interfaceBothOrderedInt
    ("BothOrderedInt",
     "If an emission should be fully ordered with both its parents."
     "otherwise it will be weighted according to how close they are.",
     &Emitter::theBothOrderedInt, false, true, false);
  static SwitchOption interfaceBothOrderedIntTrue
    (interfaceBothOrderedInt,"True","both parents fully ordered.",true);
  static SwitchOption interfaceBothOrderedIntFalse
    (interfaceBothOrderedInt,"False","weighted by distance.",false);

  static Switch<Emitter,bool> interfaceBothOrderedFS
    ("BothOrderedFS",
     "If an emission should be fully ordered with both its parents."
     "otherwise it will be weighted according to how close they are.",
     &Emitter::theBothOrderedFS, false, true, false);
  static SwitchOption interfaceBothOrderedFSTrue
    (interfaceBothOrderedFS,"True","both parents fully ordered.",true);
  static SwitchOption interfaceBothOrderedFSFalse
    (interfaceBothOrderedFS,"False","weighted by distance.",false);

  static Switch<Emitter,int> interfaceRangeMode
    ("RangeMode",
     "How the range of coherenent emissions is determined.",
     &Emitter::theRangeMode, 0, true, false);
  static SwitchOption interfaceRangeModeMin
    (interfaceRangeMode,
     "Min",
     "minimum of half mother dipole and distance to emission. Default.",
     0);
  static SwitchOption interfaceRangeModeMax
    (interfaceRangeMode,
     "Mother",
     "Half distance of mother dipole, no matter how close the emission is.",
     1);

  static Switch<Emitter,int> interfaceMinusOrderingMode
    ("MinusOrderingMode",
     "Sets how the ordering in p- is done in the virtual cascade, in relation to effective partons mainly.",
     &Emitter::theMinusOrderingMode, 0, true, false);
  static SwitchOption interfaceMinusOrderingModeEffectiveParton
    (interfaceMinusOrderingMode,
     "EffectiveParton",
     "Uses the momentum of the effective parton to order, both plus and pt.",
     0);
  static SwitchOption interfaceMinusOrderingModeEffectivePT
    (interfaceMinusOrderingMode,
     "EffectivePT",
     "Uses the pt of the effective parton, but the rapidity of the single emitting parton.",
     1);
  static SwitchOption interfaceMinusOrderingModeTrueShadow
    (interfaceMinusOrderingMode,
     "TrueShadow",
     "Use the full shadow mechanism for resolved emissions to get the "
     "incoming propagator.",
     2);
  static SwitchOption interfaceMinusOrderingModeOrderedShadow
    (interfaceMinusOrderingMode,
     "OrderedShadow",
     "Use the full shadow mechanism for resolved emissions to get the "
     "incoming propagator and checking ordering with previous emissions.",
     3);

}

