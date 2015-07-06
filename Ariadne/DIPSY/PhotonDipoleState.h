// -*- C++ -*-
#ifndef DIPSY_PhotonDipoleState_H
#define DIPSY_PhotonDipoleState_H
//
// The PhotonDipoleState class is a simple wrapper sub-class of
// DipoleState for easy construction of an initial simple q-qbar
// dipole.
//

#include "DipoleState.h"
#include "PhotonWFInfo.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the PhotonDipoleState class.
 */
class PhotonDipoleState: public DipoleState {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline PhotonDipoleState();

  /**
   * The standard constructor taking an event handler, the total
   * positive light-cone momentum and a PhotonWFInfo as argument.
   */
  PhotonDipoleState(const DipoleEventHandler & eh, Energy plus, Energy minus,
		    Ptr<PhotonWFInfo>::pointer wfi, double weight);

  /**
   * The copy constructor.
   */
  inline PhotonDipoleState(const PhotonDipoleState &);

  /**
   * The destructor.
   */
  virtual ~PhotonDipoleState();
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PhotonDipoleState & operator=(const PhotonDipoleState &);

};

}

#include "PhotonDipoleState.icc"

#endif /* DIPSY_PhotonDipoleState_H */
