// -*- C++ -*-
#ifndef Ariadne5_ColourResonanceModel_H
#define Ariadne5_ColourResonanceModel_H
//
// This is the declaration of the ColourResonanceModel class.
//

#include "ColourResonanceModel.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "Ariadne/Cascade/ResonanceParton.h"
#include "Ariadne/Cascade/Emission.h"
#include "PseudoParton.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * ColourResonanceModel is a helper class to be used when emitting
 * from a coloured resonance.
 *
 * @see \ref ColourResonanceModelInterfaces "The interfaces"
 * defined for ColourResonanceModel.
 */
class ColourResonanceModel: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ColourResonanceModel();

  /**
   * The destructor.
   */
  virtual ~ColourResonanceModel();
  //@}

public:

  /**
   * Check if the given resonance product should still be considered
   * special, or if it can be treated as any other parton.
   */
  virtual bool stillSpecial(tResParPtr p) const;

  /**
   * Return a pseudo particle to be used in normal dipole emission
   * from the given resonance product.
   */
  virtual PseudoParton getPseudoParton(tResParPtr p) const;

  /**
   * If the given resonance product is still special, generate
   * possible gluon emissions from the internal coloue line. If \a
   * iside is true thie is for the iPart() of the emitting dipole.
   */
  virtual EmPtr
  generateInternalGluon(const EmitterBase &, const QCDDipole & dip,
			tResParPtr p) const;

  /**
   * Set momentum after emission from a pseudo parton.
   */
  virtual void setMomentum(PseudoParton pp, const Lorentz5Momentum & p) const;

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ColourResonanceModel & operator=(const ColourResonanceModel &);

};

}

#endif /* Ariadne5_ColourResonanceModel_H */
