// -*- C++ -*-
#ifndef DIPSY_VirtualPhoton_H
#define DIPSY_VirtualPhoton_H
//
// This is the declaration of the VirtualPhoton class.
//

#include "Ariadne/DIPSY/VectorMesonBase.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The VirtualPhoton class represents the perturbative
 * quark--anti-quark dipole wave function of a virtual photon with a
 * given virtuality.
 *
 * @see \ref VirtualPhotonInterfaces "The interfaces"
 * defined for VirtualPhoton.
 */
class VirtualPhoton: public VectorMesonBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline VirtualPhoton()
    : theQ2(0.0*GeV2), theShrinkR(0.0*InvGeV),thePolarisation(0),
      theVMDLightA(0.0), theVMDLightR(0.0*InvGeV),
      theVMDLightW(0.0*InvGeV), theVMDCharmA(0.0), theVMDCharmR(0.0*InvGeV),
      theVMDCharmW(0.0*InvGeV), rMax(ZERO), r2Psi2Int(ZERO), r2Psi2Max(ZERO) {}

  /**
   * The copy constructor.
   */
  inline VirtualPhoton(const VirtualPhoton & x)
    : VectorMesonBase(x), theQ2(x.theQ2), theShrinkR(x.theShrinkR),
      thePolarisation(x.thePolarisation), theVMDLightA(x.theVMDLightA),
      theVMDLightR(x.theVMDLightR), theVMDLightW(x.theVMDLightW),
      theVMDCharmA(x.theVMDCharmA), theVMDCharmR(x.theVMDCharmR),
      theVMDCharmW(x.theVMDCharmW), rMax(x.rMax),
      r2Psi2Int(x.r2Psi2Int), r2Psi2Max(x.r2Psi2Max) {}

  /**
   * The destructor.
   */
  virtual ~VirtualPhoton();
  //@}

public:

  /** @name Main virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Initialize the wavefunction for the given DipoleEventHandler.
   */
  virtual void initialize(const DipoleEventHandler &);

  /**
   * Generate a dipole state according to this wavefunction, given an
   * event handler and a positive light-cone momentum component.
   */
  virtual DipoleStatePtr generate(const DipoleEventHandler & eh, Energy plus);

  /**
   * Generate a dipole state according to this wavefunction, given an
   * event handler and a positive light-cone momentum component.
   */
  virtual DipoleStatePtr generate2(const DipoleEventHandler & eh, Energy plus);

  /**
   * Return (the absolute value of) the wave function for the given
   * polarization \a pol, quark and antiquark helicities, \a h and \a
   * hbar, quark flavour \a flav, energy fraction \a z and transverse
   * size \a r.
   */
  virtual Energy psi(int pol, int h, int hbar, int flav,
		     InvEnergy r, double z);

  /**
   * The purely perturbative version of the wavefunction.
   */
  virtual Energy pertPsi(int pol, int h, int hbar, int flav,
			 InvEnergy r, double z);

  /**
   * r*2*pi*(sum_i Psi_i^2(r,z)) with sum i over helicities. this describes
   * the weight to get a dipole with certain r and z.
   **/
  virtual Energy sumPsi2(InvEnergy r, double z);

  /**
   * returns the correction on the wavefunction due to vector meson resonance
   * at a given dipole size for a given quark flavour.
   **/
  virtual double VMDCorr(InvEnergy r, int flav);

  /**
   * Return the invariant mass squared of the particle.
   */
  virtual Energy2 m2() const;
  //@}

  /**
   * The virtuality of the photon.
   */
  inline Energy2 Q2() const {
    return theQ2;
  }

  /**
   * The confinement range.
   */
  inline InvEnergy shrinkR() const {
    return theShrinkR;
  };

  /**
   * The amplitude of the vector meson resonance for light quarks.
   */
  inline double VMDLightA() const {
    return theVMDLightA;
  };

  /**
   * A flag showing which polarisation the photon has.
   */
  inline int polarisation() const {
    return thePolarisation;
  };

  /**
   * The typical dipole size of the vector meson resonance for light quarks.
   */
  inline InvEnergy VMDLightR() const {
    return theVMDLightR;
  };

  /**
   * The width in dipole size of the vector meson resonance for light quarks.
   */
  inline InvEnergy VMDLightW() const {
    return theVMDLightW;
  };

  /**
   * The amplitude of the vector meson resonance for charm quarks.
   */
  inline double VMDCharmA() const {
    return theVMDCharmA;
  };

  /**
   * The typical dipole size of the vector meson resonance for charm quarks.
   */
  inline InvEnergy VMDCharmR() const {
    return theVMDCharmR;
  };

  /**
   * The width in dipole size of the vector meson resonance for charm quarks.
   */
  inline InvEnergy VMDCharmW() const {
    return theVMDCharmW;
  };

  /**
   * Set the corresponding particle. Check that it is actually a
   * photon.
   */
  virtual void setParticle(PDPtr);

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const {
    return new_ptr(*this);
  }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {
    return new_ptr(*this);
  }
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The virtuality of the photon.
   */
  Energy2 theQ2;

  /**
   * Parameter for the confinement range
   */
  InvEnergy theShrinkR;

  /**
   * 0: both T and L, 1: Transverse only, 2: Longitudinal only
   */
  int thePolarisation;

  /**
   * Parameter for the vector meson resonance for light and charm quarks.
   */
  double theVMDLightA;
  InvEnergy theVMDLightR;
  InvEnergy theVMDLightW;
  double theVMDCharmA;
  InvEnergy theVMDCharmR;
  InvEnergy theVMDCharmW;

  /**
   * The maximum dipole size that will be considered.
   */
  InvEnergy rMax;

  /**
   * The integral of r^2 times the sum of the square of all wavefunctions, up to rMax.
   **/
  InvEnergy2 r2Psi2Int;

  /**
   * The max of r^2 times the sum of the square of all the wavefunctions, up to rMax.
   **/
  InvEnergy r2Psi2Max;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VirtualPhoton & operator=(const VirtualPhoton &);

};

}

#endif /* DIPSY_VirtualPhoton_H */
