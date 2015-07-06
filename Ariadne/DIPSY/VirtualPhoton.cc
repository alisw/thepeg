// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VirtualPhoton class.
//

#include "VirtualPhoton.h"
#include "PhotonWFInfo.h"
#include "ThePEG/Interface/Switch.h"
#include "PhotonDipoleState.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/PDT/EnumParticles.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "gsl/gsl_sf_bessel.h"

using namespace DIPSY;

VirtualPhoton::~VirtualPhoton() {}

void VirtualPhoton::initialize(const DipoleEventHandler &) {

  cout << "preparing virtual photon generation for " << this << endl;
  //find the normalisation and maximum value of r^3*2*pi*(sum_i Psi_i^2(r,z))
  //this is the distribution that will be used to generate dipoles, and it is
  //only the wavefunction, times an extra r^2 for interaction probability
  //and converging integral
  rMax = 10.0/GeV;
  r2Psi2Int = ZERO;
  r2Psi2Max = ZERO;
  int samples = 10000;
  InvEnergy sum = ZERO;
  for ( int i = 0; i < samples; i++ ) {
    double z = UseRandom::rnd();
    InvEnergy r = UseRandom::rnd()*rMax;
    InvEnergy value = sqr(r)*sumPsi2(r, z);
    if ( value > r2Psi2Max ) r2Psi2Max = value;
    sum += value;
  }
  r2Psi2Int = sum*rMax/double(samples);
  cout << "photon generation ready. Max value: " << r2Psi2Max*GeV << "GeV-1, integral: "
       << r2Psi2Int*GeV2 << "GeV-2" << endl;

}

Energy2 VirtualPhoton::m2() const {
  return -Q2();
}

DipoleStatePtr VirtualPhoton::
generate2(const DipoleEventHandler & eh, Energy plus) {

  //first generate an r and z with (normalised) distribution P_0(r) = r^2 P_true(r)
  InvEnergy r;
  double z;
  while (true) {
    r = UseRandom::rnd()*rMax;
    z = UseRandom::rnd();
    if ( sqr(r)*sumPsi2(r, z) > UseRandom::rnd()*r2Psi2Max ) break;
  }

  //set weight to account for normalisation of P_0(r), and to make up
  //for the extra r^2 to get to P_true(r)
  double weight = r2Psi2Int/sqr(r);

  //pick flavour according to the square of the wavefunction for this r and z
  Selector<int,Energy2> selFlav;
  for ( int f = 1; f < maxFlav(); ++f ) {
    Energy2 prob = ZERO;
    if ( polarisation() == 0 || polarisation() == 2 )
      prob += sqr(psi(0, 1, 1, f, r, z))*2.0;
    if ( polarisation() == 0 || polarisation() == 1 ) {
      prob += sqr(psi(1, 1, -1, f, r, z));
      prob += sqr(psi(1, -1, 1, f, r, z));
      prob += sqr(psi(1, 1, 1, f, r, z));
    }
    selFlav.insert((f%2? 1.0: 4.0)*prob, f);
  }
  int f = selFlav[UseRandom::rnd()];

  // Now pick helicity, with the weights of the squared wavefunction
  //for this r, z and f.
  Energy2 PLhh = sqr(psi(0, 1, 1, f, r, z))*2.0;
  Energy2 PTpm = sqr(psi(1, 1, -1, f, r, z));
  Energy2 PTmp = sqr(psi(1, -1, 1, f, r, z));
  Energy2 PThh = sqr(psi(1, 1, 1, f, r, z));
  if ( polarisation() == 1 )
    PLhh = ZERO;
  if ( polarisation() == 2 ) {
    PTpm = ZERO;
    PTmp = ZERO;
    PThh = ZERO;
  }
  Energy2 sum = PLhh + PTpm + PTmp + PThh;
  Energy2 rsum = sum*UseRandom::rnd();
  Ptr<PhotonWFInfo>::pointer wfi;
  if ( PLhh > rsum ) {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 0, 1, 1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, 0, -1, -1, f));
  }
  else if ( PLhh + PThh > rsum ) {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 1, 1, 1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, -1, -1, -1, f));
  }
  else if ( PLhh + PThh + PTpm > rsum ) {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 1, 1, -1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, -1, 1, -1, f));
  }
  else {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 1, -1, 1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, -1, -1, 1, f));
  }
  
  return new_ptr(PhotonDipoleState(eh, plus, -Q2()/plus, wfi, weight));


}

DipoleStatePtr VirtualPhoton::
generate(const DipoleEventHandler & eh, Energy plus) {
  // z is always chosen uniformly
  double z = UseRandom::rnd();

  // flavour is chosen with (normalized) integrated probability
  // proportional to ef^2 exp(-2*epsf*r)
  Selector<int,InvEnergy> sel;
  Energy eps = sqrt(z*(1.0 - z)*Q2() + sqr(qmass[abs(maxFlav())]));
  sel.insert((abs(maxFlav())%2? 1.0: 4.0)/(18.0*eps), abs(maxFlav()));
  for ( int f = 1; f < maxFlav(); ++f ) {
    Energy eps = sqrt(z*(1.0 - z)*Q2() + sqr(qmass[f]));
    sel.insert((f%2? 1.0: 4.0)/(18.0*eps), f);
  }
  int f = sel[UseRandom::rnd()];

  // r is chosen according to (normalized) probability exp(-2*epsf*r).
  eps = sqrt(z*(1.0 - z)*Q2() + sqr(qmass[f]));
  InvEnergy r =-log(UseRandom::rnd())/(2.0*eps);

  // Przf is the probability with which we have chosen r, z and f
  Energy Przf = (f%2? 1.0: 4.0)*exp(-2.0*r*eps)/(9.0*sel.sum());

  // Now find the squared wavefunction (differential in r) for the
  // different helicities and choose one.
  Energy PLhh = r*sqr(psi(0, 1, 1, f, r, z))*4.0*Constants::pi;
  Energy PTpm = r*sqr(psi(1, 1, -1, f, r, z))*2.0*Constants::pi;
  Energy PTmp = r*sqr(psi(1, -1, 1, f, r, z))*2.0*Constants::pi;
  Energy PThh = r*sqr(psi(1, 1, 1, f, r, z))*2.0*Constants::pi;
  if ( polarisation() == 1 )
    PLhh = ZERO;
  if ( polarisation() == 2 ) {
    PTpm = ZERO;
    PTmp = ZERO;
    PThh = ZERO;
  }
  Energy sum = PLhh + PTpm + PTmp + PThh;
  Energy rsum = sum*UseRandom::rnd();
  Ptr<PhotonWFInfo>::pointer wfi;
  double weight = sum/Przf;
  if ( PLhh > rsum ) {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 0, 1, 1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, 0, -1, -1, f));
  }
  else if ( PLhh + PThh > rsum ) {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 1, 1, 1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, -1, -1, -1, f));
  }
  else if ( PLhh + PThh + PTpm > rsum ) {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 1, 1, -1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, -1, 1, -1, f));
  }
  else {
    if ( UseRandom::rndbool() )
      wfi = new_ptr(PhotonWFInfo(this, r, z, 1, -1, 1, f));
    else
      wfi = new_ptr(PhotonWFInfo(this, r, z, -1, -1, 1, f));
  }

  // cout << setprecision(7);
  // if ( r == ZERO ) cout << "r is identically zero" << endl;
  // cout << "generated photon with r = " << r*GeV << ", z = " << z
  //      << ", weight = " << weight << endl;
  
  return new_ptr(PhotonDipoleState(eh, plus, -Q2()/plus, wfi, weight));
}

Energy VirtualPhoton::
pertPsi(int pol, int h, int hbar, int flav, InvEnergy r, double z) {
  double ef = flav%2? -1.0/3.0: 2.0/3.0;
  Energy eps = sqrt(z*(1.0 - z)*Q2() + sqr(qmass[flav]));
  double epsr = eps*r;
  if ( epsr > 100.0 ) return ZERO;
  if ( epsr == 0.0 ) return ZERO;
  if ( pol == 0 ) {
    if ( h != hbar ) return 0.0*GeV;
    return sqrt(SM().alphaEM()*SM().Nc()*Q2())*z*(1.0 - z)*ef*
      gsl_sf_bessel_K0(epsr)/Constants::pi;
  } else {
    if ( h == pol && hbar == -pol )
      return sqrt(SM().alphaEM()*SM().Nc()/2.0)*ef*z*eps*
	gsl_sf_bessel_K1(epsr)/Constants::pi;
    else if ( h == -pol && hbar == pol )
      return sqrt(SM().alphaEM()*SM().Nc()/2.0)*ef*(1.0 - z)*eps*
	gsl_sf_bessel_K1(epsr)/Constants::pi;
    else if ( h == pol && hbar == pol )
      return sqrt(SM().alphaEM()*SM().Nc()/2.0)*ef*qmass[flav]*
	gsl_sf_bessel_K0(epsr)/Constants::pi;
    else
      return 0.0*GeV;
  }
}

// int status =  gsl_sf_bessel_K1_e(epsr,&result);
// if (status) {
//   //handle error.
//  }
// result.val;
// typedef struct
// {
//   double val;
//   double err;
// } gsl_sf_result;

Energy VirtualPhoton::
psi(int pol, int h, int hbar, int flav, InvEnergy r, double z) {
  // Start out with the perturbative wavefunciton for the unconfined r.
  InvEnergy rPert = r;
  double jacobian = 1.0;
  //taking confinement into account, what perturbative (nonconfined)
  //r_pert is the final confined r_conf associated with?
  if ( shrinkR() != ZERO ) rPert = shrinkR()*sqrt(-Math::exp1m(sqr(r/shrinkR())));
  //events are generated in dr_conf, but the wavefuntion create them
  //with dr_pert, so we need an extra dr_pert/dr_conf to account for how
  //the dip sizes gets clumped up around confinement scale.
  if ( shrinkR() != ZERO ) 
    jacobian = r/shrinkR()*exp(sqr(r/shrinkR()))/sqrt(expm1(sqr(r/shrinkR())));

  Energy ret = pertPsi(pol, h, hbar, flav, rPert, z);
  //vector meson resonance enhances the wavefunction if the final confined r
  //is about the size of a vector emson.
  ret *= sqrt(VMDCorr(r, flav)*jacobian);

  return ret;
}

Energy VirtualPhoton::
sumPsi2(InvEnergy r, double z) {
  //sum all wavefunctions contribution to this r and z
  Energy ret = ZERO;
  //sum over flavours
  for ( int f = 1; f < maxFlav(); ++f ) {
    //sum over helicities
    if ( polarisation() == 0 || polarisation() == 2 )
      ret += r*sqr(psi(0, 1, 1, f, r, z))*4.0*Constants::pi;
    if ( polarisation() == 0 || polarisation() == 1 ) {
      ret += r*sqr(psi(1, 1, -1, f, r, z))*2.0*Constants::pi;
      ret += r*sqr(psi(1, -1, 1, f, r, z))*2.0*Constants::pi;
      ret += r*sqr(psi(1, 1, 1, f, r, z))*2.0*Constants::pi;
    }
  }
  return ret;
}

double VirtualPhoton::VMDCorr(InvEnergy r, int flav) {
  if ( VMDLightA() == 0.0 ) return 1.0;
  if ( flav == 1 || flav == 2 || flav == 3 )
    return (1.0 + VMDLightA()*exp(-sqr((r-VMDLightR())/VMDLightW())))/
      (1.0 + VMDLightA()*exp(-sqr(VMDLightR()/VMDLightW())));
  else if ( flav == 4 )
    return (1.0 + VMDCharmA()*exp(-sqr((r-VMDCharmR())/VMDCharmW())))/
      (1.0 + VMDCharmA()*exp(-sqr(VMDCharmR()/VMDCharmW())));

  cerr << "vector meson dominance parameters only added for up, down, strange and charm. Heavier quarks get the charm parameters. flavour is " << flav << endl;
  return (1.0 + VMDCharmA()*exp(-sqr((r-VMDCharmR())/VMDCharmW())))/
    (1.0 + VMDCharmA()*exp(-sqr(VMDCharmR()/VMDCharmW())));

}

void VirtualPhoton::setParticle(PDPtr p) {
  if ( p->id() != ParticleID::gamma )
    throw InterfaceException()
      << "Cannot set " << p->name()
      << " as particle for a VirtualPhoton wave function. "
      << "Only Photons are allowed." << Exception::warning;
  WaveFunction::setParticle(p);
}

void VirtualPhoton::persistentOutput(PersistentOStream & os) const {
  os << ounit(theQ2, GeV2) << ounit(theShrinkR, InvGeV) << thePolarisation << theVMDLightA
     << ounit(theVMDLightR, InvGeV) << ounit(theVMDLightW, InvGeV) << theVMDCharmA
     << ounit(theVMDCharmR, InvGeV) << ounit(theVMDCharmW, InvGeV) << ounit(rMax, InvGeV)
     << ounit(r2Psi2Int, sqr(InvGeV)) << ounit(r2Psi2Max, InvGeV);
}

void VirtualPhoton::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theQ2, GeV2) >> iunit(theShrinkR, InvGeV) >> thePolarisation >> theVMDLightA
     >> iunit(theVMDLightR, InvGeV) >> iunit(theVMDLightW, InvGeV) >> theVMDCharmA
     >> iunit(theVMDCharmR, InvGeV) >> iunit(theVMDCharmW, InvGeV) >> iunit(rMax, InvGeV)
     >> iunit(r2Psi2Int, sqr(InvGeV)) >> iunit(r2Psi2Max, InvGeV);
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<VirtualPhoton,DIPSY::VectorMesonBase>
  describeDIPSYVirtualPhoton("DIPSY::VirtualPhoton", "libAriadne5.so libDIPSY.so");


void VirtualPhoton::Init() {

  static ClassDocumentation<VirtualPhoton> documentation
    ("The VirtualPhoton class represents the perturbative quark--anti-quark "
     "dipole wave function of a virtual photon with a given virtuality.");

  static Parameter<VirtualPhoton,Energy2> interfaceQ2
    ("Q2",
     "The virtuality of the photon in units of GeV\\f$^2\\f$.",
     &VirtualPhoton::theQ2, GeV2, 0.0*GeV2, 0.0*GeV2, 0*GeV2,
     true, false, Interface::lowerlim);

  static Parameter<VirtualPhoton, InvEnergy> interfaceShrinkR
    ("ShrinkR",
     "Sets the confinement scale of the photon wavefunction. Due to different"
"parametrisations, this paramter should be about 1.5 times rmax in DipoleEventHandler.",
     &VirtualPhoton::theShrinkR, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<VirtualPhoton, double> interfaceVMDLightA
    ("VMDLightA",
     "Sets the amplitude of the vector meson resonence boost for the wavefunction.",
     &VirtualPhoton::theVMDLightA, 0.0, 0.0, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Parameter<VirtualPhoton, InvEnergy> interfaceVMDLightR
    ("VMDLightR",
     "Sets the typical size of the vector meson resonence in wavefunction.",
     &VirtualPhoton::theVMDLightR, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<VirtualPhoton, InvEnergy> interfaceVMDLightW
    ("VMDLightW",
     "Sets the width of the vector meson resonence in wavefunction.",
     &VirtualPhoton::theVMDLightW, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<VirtualPhoton, double> interfaceVMDCharmA
    ("VMDCharmA",
     "Sets the amplitude of the vector meson resonence boost for the wavefunction.",
     &VirtualPhoton::theVMDCharmA, 0.0, 0.0, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Parameter<VirtualPhoton, InvEnergy> interfaceVMDCharmR
    ("VMDCharmR",
     "Sets the typical size of the vector meson resonence in wavefunction.",
     &VirtualPhoton::theVMDCharmR, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<VirtualPhoton, InvEnergy> interfaceVMDCharmW
    ("VMDCharmW",
     "Sets the width of the vector meson resonence in wavefunction.",
     &VirtualPhoton::theVMDCharmW, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Switch<VirtualPhoton,int> interfacePolarisation
    ("Polarisation",
     "What helicity of the photon to run.",
     &VirtualPhoton::thePolarisation, 0, true, false);
  static SwitchOption interfaceHelicityBoth
    (interfacePolarisation,
     "Both",
     "Use both T and L helicities, effectiviely giving F_2.",
     0);
  static SwitchOption interfaceHelicityTransverse
    (interfacePolarisation,
     "Transverse",
     "Use only transverse helicity photons.",
     1);
  static SwitchOption interfaceHelicityLongitudinal
    (interfacePolarisation,
     "Longitudinal",
     "Use only longitudinal polarisation, giving F_L",
     2);
}

