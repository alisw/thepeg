// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BudnevPDF class.
//

#include "BudnevPDF.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;

BudnevPDF::BudnevPDF() 
  : _q2min(ZERO), _q2max(2.*GeV2),_q02(0.71*GeV2),_mup2(7.78) 
{}

bool BudnevPDF::canHandleParticle(tcPDPtr particle) const {
  return ( abs(particle->id()) == ParticleID::pplus );
}

cPDVector BudnevPDF::partons(tcPDPtr) const {
  // only photon
  return cPDVector(1,getParticleData(ParticleID::gamma));
}

double BudnevPDF::xfl(tcPDPtr proton, tcPDPtr gamma, Energy2 qq,
                      double l, Energy2 ) const {

  if(gamma->id()!=ParticleID::gamma) return 0.;
  double x(exp(-l));
  
  // photon virtuality allowed by kinematics
  Energy2 qqkinmin = sqr(proton->mass()*x)/(1-x);
  
  //electric form factor for proton
  double fe = (4*sqr(proton->mass())*ge2(qq) + qq*gm2(qq))/(4*sqr(proton->mass()) + qq);
  //magnetic form factor
  double fm = gm2(qq);
  
  return SM().alphaEM()/Constants::pi*((1.-x)*(1-qqkinmin/qq)*fe + 0.5*sqr(x)*fm);
}

double BudnevPDF::xfvl(tcPDPtr, tcPDPtr, Energy2, double,
				   Energy2) const {
  // valence density is zero
  return 0.0;
}



// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BudnevPDF,PDFBase>
describeBudnevPDF("ThePEG::BudnevPDF", "BudnevPDF.so");

void BudnevPDF::Init() {

  static ClassDocumentation<BudnevPDF> documentation
    ("The BudnevPDF provides the PDF for a photon inside"
     " an incoming photon using the results of Phys. Rept. 15 (1975) 181.",
     "The BudnevPDF provides the PDF for a photon inside"
     " an incoming photon based on the results of \\cite{Budnev:1974de} was used.",
     "\\bibitem{Budnev:1974de}"
     "V.~M.~Budnev, I.~F.~Ginzburg, G.~V.~Meledin and V.~G.~Serbo,"
     "%``The Two photon particle production mechanism. Physical problems. Applications. Equivalent photon approximation,''"
     "Phys.\\ Rept.\\  {\\bf 15} (1975) 181."
     "doi:10.1016/0370-1573(75)90009-5"
     "%%CITATION = doi:10.1016/0370-1573(75)90009-5;%%"
     "%1212 citations counted in INSPIRE as of 27 Nov 2019");
     
  static Parameter<BudnevPDF,Energy2> interfaceQ2Min
    ("Q2Min",
     "Minimum value of the magnitude of Q^2 for the photon",
     &BudnevPDF::_q2min, GeV2, ZERO, ZERO, 100.0*GeV2,
     false, false, Interface::limited);

  static Parameter<BudnevPDF,Energy2> interfaceQ2Max
    ("Q2Max",
     "Maximum value of the magnitude of Q^2 for the photon",
     &BudnevPDF::_q2max, GeV2, 4.0*GeV2, ZERO, 100.0*GeV2,
     false, false, Interface::limited);

}

double BudnevPDF::
flattenScale(tcPDPtr proton, tcPDPtr, const PDFCuts & c,
	     double l, double z, double & jacobian) const {
  double x = exp(-l);
  
  Energy2 qqmax = min(_q2max,0.25*sqr(x)*c.sMax());
  Energy2 qqmin = max(_q2min, sqr(proton->mass()*x)/(1-x));
  if(qqmin>=qqmax) {
    jacobian = 0.;
    return 0.;
  }
  double low(log(qqmin/c.scaleMaxL(l))),upp(log(qqmax/c.scaleMaxL(l)));
  // jacobian factor
  jacobian *= log(qqmax/qqmin);
  return exp(low+z*(upp-low));
}

double BudnevPDF::flattenL(tcPDPtr, tcPDPtr, const PDFCuts & c,
			 double z, double & jacobian) const {
  jacobian *= c.lMax() - c.lMin();
  return c.lMin() + z*(c.lMax() - c.lMin());
}

void BudnevPDF::persistentOutput(PersistentOStream & os) const {
  os << ounit(_q2min,GeV2) << ounit(_q2max,GeV2);
}

void BudnevPDF::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_q2min,GeV2) >> iunit(_q2max,GeV2);
}
