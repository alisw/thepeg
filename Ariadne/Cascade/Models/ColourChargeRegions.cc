// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourChargeRegions class.
//

#include "ColourChargeRegions.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Current.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;


ColourChargeRegions::ColourChargeRegions() {}

ColourChargeRegions::~ColourChargeRegions() {}

IBPtr ColourChargeRegions::clone() const {
  return new_ptr(*this);
}

IBPtr ColourChargeRegions::fullclone() const {
  return new_ptr(*this);
}

double ColourChargeRegions::preweight(const Emission & e) const {
  tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(e.cdipole);
  if ( !d || d->iPart()->isG() || d->oPart()->isG() ) return 1.0;
  return 9.0/8.0;
}

double ColourChargeRegions::reweight(const Emission & e) const {
  tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(e.cdipole);
  if ( !d ) return 1.0;
  if ( &e != lastEmission ) {
    lastEmission = &e;
    setup(*d);
  }
  return charge(e.y)/(3.0/2.0);
}

void ColourChargeRegions::setup(const QCDDipole & d) const {
  static const double Cg = 3.0/2.0;
  static const double Cq = 4.0/3.0;

  colregions.clear();
  acoregions.clear();

  Energy rhomin = Current<AriadneHandler>()->pTCut();
  double ymax = log(sqrt(d.sdip())/rhomin);

  // First check anti-colour side for rapidity regions.
  tcParPtr op = d.oPart();
  double ymaxl = ymax; // The current max rapidity
  double yl = ymaxl; // The last rapidity where we have inserted a fold. 
  // Set the limit - close to the parton we know the colour charge
  acoregions[yl] =  op->isG()? Cg: Cq;
  while ( op->emission() ) {
    // Now check where this parton was emitted. Clearly only the fold
    // this parton produced will carry its charge.
    double y = ymaxl - log(op->emission()->rho/rhomin);
    // We now need to increase the maximum rapidity as we go to the
    // parent. This is the distance to the parent minus the extra fold
    // from the emission.
    ymaxl = y + (op->emission()->ymax + op->emission()->y);
    op = op->emission()->antiColourParent;
    if ( y < yl ) {
      // Only insert a new fold limit if it is ordered
      acoregions[y] = op->isG()? Cg: Cq;
      yl = y;
    } else {
      acoregions.begin()->second = op->isG()? Cg: Cq;
    }
  }
  // *** ATTENTION **** The following is maybe wrong. In fact the
  // whole thing may need to be reworked.
  if ( op->orig() && op->orig()->scale() > sqr(rhomin) ){
    double y = ymaxl - log(sqrt(op->orig()->scale())/rhomin); 
    if ( y < yl ) acoregions[y] = 0.0;
  }

  // Then check the colour side.
  tcParPtr ip = d.iPart();
  ymaxl = ymax;
  yl = ymax;
  colregions[yl] = ip->isG()? Cg: Cq;
  while ( ip->emission() ) {
    double y = ymaxl - log(ip->emission()->rho/rhomin);
    ymaxl = y + (ip->emission()->ymax - ip->emission()->y);
    ip = ip->emission()->colourParent;
    if ( y < yl ) {
      colregions[y] = ip->isG()? Cg: Cq;
      yl = y;
    } else {
      colregions.begin()->second = ip->isG()? Cg: Cq;
    }
  }
  if ( ip->orig() && ip->orig()->scale() > sqr(rhomin) ) {
    double y = ymaxl - log(sqrt(ip->orig()->scale())/rhomin); 
    if ( y < yl ) colregions[y] = 0.0;
  }
}

double ColourChargeRegions::charge(double y) const {
  map<double,double>::const_iterator colit = colregions.upper_bound(-y);
  map<double,double>::const_iterator acoit = acoregions.upper_bound(y);


  // If we are outside the reach of either original parton let the
  // other decide the charge. Or if outside the reach of both, take
  // the geometric mean.
  if ( colit->second == 0.0 && acoit->second == 0.0 )
    return sqrt((++colit)->second * (++acoit)->second);
  if ( colit->second == 0.0 ) return acoit->second;
  if ( acoit->second == 0.0 ) return colit->second;

  // If we are in the charge regions of non-original partons on both
  // sides, our reconstruction basicalle failed. Return the geometric
  // mean.
  if ( colit != colregions.begin() && acoit != acoregions.begin() )
    return sqrt(colit->second * acoit->second);

  // If we are in the charge region of non-original partons on one
  // side, return that charge.
  if ( colit != colregions.begin() ) return colit->second;
  if ( acoit != acoregions.begin() ) return acoit->second;

  // If we are properly inbetween original partons, interpolate with
  // rapidity.
  return ((y + colit->first)*acoit->second +
	  (y + acoit->first)*colit->second)/
   (colit->first + acoit->first);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ColourChargeRegions::persistentOutput(PersistentOStream &) const {}

void ColourChargeRegions::persistentInput(PersistentIStream &, int) {}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<ColourChargeRegions,Ariadne5::ReweightBase>
describeAriadne5ColourChargeRegions("Ariadne5::ColourChargeRegions",
			    "libAriadne5.so.so");

void ColourChargeRegions::Init() {

  static ClassDocumentation<ColourChargeRegions> documentation
    ("The ColourChargeRegions class can be used to reweight gluon "
     "emissions so that a more correct colour charge is assigned to "
     "an emission probability, based not only on the charge of the"
     " emitting partons in the dipole, but also on where in rapidity "
     "the gluon is emitted. This uses a simplified version of the "
     "model by Eden and Gustafson [hep-ph/9805228].",

     "Colour-factor corrections in gluon emissions given by a simplified "
     "version of the model in \\cite{Eden:1998ig}.",

     "\\bibitem{Eden:1998ig}\n  P.~Ed\\'en and G.~Gustafson,\n"
     "%``Energy and virtuality scale dependence in quark and gluon jets,''\n"
     "JHEP {\\bf 9809}, 015 (1998) [arXiv:hep-ph/9805228].\n"
     "%%CITATION = JHEPA,9809,015;%%\n");

}

void ColourChargeRegions::debugme() const {
  ReweightBase::debugme();
  cerr << endl << fullName() << endl
       << "** Outgoing colour line: " << setprecision(3) << endl;
  for ( map<double,double>::const_reverse_iterator it = colregions.rbegin();
	it != colregions.rend(); ++it )
    cerr << "   " << it->first << " " << it->second << endl;
  cerr << "** Incoming anti-colour line (in reverse): " << endl;
  for ( map<double,double>::const_iterator it = acoregions.begin();
	it != acoregions.end(); ++it )
    cerr << "   " << -it->first << " " << it->second << endl;
  if ( !lastEmission ) return;
  tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(lastEmission->cdipole);
  if ( !d ) return;
  Energy rhomin = Current<AriadneHandler>()->pTCut();
  double ymax = log(sqrt(d->sdip())/rhomin);
  cerr << "** This is how we got here:" << endl;
  cerr << "   Maximum rapidity: " << ymax << endl;
  cerr << " * Colour parton:" << endl;
  tcParPtr ip = d->iPart();
  ip->debug();
  while ( ip->emission() ) {
    cerr << "   emitted at kappa/2 = " << log(ip->emission()->rho/rhomin)
	 << " and y = " << ip->emission()->y << " (max = "
	 << ip->emission()->ymax << ") from" << endl;
    ip = ip->emission()->colourParent;
    ip->debug();
  }
  cerr << " * Anti-colour parton:" << endl;
  tcParPtr op = d->oPart();
  op->debug();
  while ( op->emission() ) {
    cerr << "   emitted at kappa/2 = " << log(op->emission()->rho/rhomin)
	 << " and y = " << op->emission()->y << " (max = "
	 << op->emission()->ymax << ") from" << endl;
    op = op->emission()->antiColourParent;
    op->debug();
  }
}
