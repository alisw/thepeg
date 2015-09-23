// -*- C++ -*-
#ifndef DIPSY_PhotonWFInfo_H
#define DIPSY_PhotonWFInfo_H
//
// This is the declaration of the PhotonWFInfo class.
//

#include "Ariadne/DIPSY/WFInfo.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The PhotonWFInfo class describes the polarizations and energy
 * sharing of a DipoleState constructed from a VirtualPhoton
 * wavefunction.
 */
class PhotonWFInfo: public WFInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor, giving the size, \a rini, energy
   * sharing, \a zini, the photon polarization, \a polini, the quark
   * and anti-quark helicities, \a hini and \a hbarini, and the quark
   * flavour, \a flini, as arguments
   */
  inline PhotonWFInfo(tcWaveFunctionPtr wfin = tcWaveFunctionPtr(),
		      InvEnergy rini = 0.0*InvGeV, double zini = 0.0,
		      int polini = 0, int hini = 0, int hbarini = 0,
		      int flini = 0)
    : WFInfo(wfin, rini), theZ(zini), thePol(polini),
      theH(hini), theHbar(hbarini), theFlav(flini) {}

  /**
   * The copy constructor.
   */
  inline PhotonWFInfo(const PhotonWFInfo & x)
    : WFInfo(x), theZ(x.theZ), thePol(x.thePol),
      theH(x.theH), theHbar(x.theHbar), theFlav(x.theFlav) {}
  
  /**
   * The destructor.
   */
  virtual ~PhotonWFInfo();
  //@}

public:

  /**
   * The energy fraction of the quark.
   */
  inline double z() const {
    return theZ;
  }

  /**
   * The polarization (-1, 0 or 1) of the photon.
   */
  inline int pol() const {
    return thePol;
  }

  /**
   * The helicity (+-1) of the quark.
   */
  inline int h() const {
    return theH;
  }

  /**
   * The helicity (+-1) of the antiquark.
   */
  inline int hbar() const {
    return theHbar;
  }

  /**
   * Get the flavour of the quark.
   */
  inline int flav() const {
    return theFlav;
  }

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

private:

  /**
   * The energy fraction of the quark.
   */
  double theZ;

  /**
   * The polarization (-1, 0 or 1) of the photon.
   */
  int thePol;

  /**
   * The helicity (+-1) of the quark.
   */
  int theH;

  /**
   * The helicity (+-1) of the antiquark.
   */
  int theHbar;

  /**
   * The flavour of the quark.
   */
  int theFlav;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PhotonWFInfo & operator=(const PhotonWFInfo &);

};

}

#endif /* DIPSY_PhotonWFInfo_H */
