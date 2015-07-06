// -*- C++ -*-
#ifndef DIPSY_VectorMesonBase_H
#define DIPSY_VectorMesonBase_H
//
// This is the declaration of the VectorMesonBase class.
//

#include "Ariadne/DIPSY/WaveFunction.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The VectorMesonBase class inherits from WaveFunction and is the
 * base class of all vector meson wave functions. It includes abstract
 * functions for different polarization and helicity components.
 *
 * @see \ref VectorMesonBaseInterfaces "The interfaces"
 * defined for VectorMesonBase.
 */
class VectorMesonBase: public WaveFunction {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline VectorMesonBase();

  /**
   * The copy constructor.
   */
  inline VectorMesonBase(const VectorMesonBase &);

  /**
   * The destructor.
   */
  virtual ~VectorMesonBase();
  //@}

public:

  /** @name Abstract functions for different wave function components. */
  //@{
  /**
   * Return (the absolute value of) the wave function for the given
   * polarization \a pol, quark and antiquark helicities, \a h and \a
   * hbar, quark flavour \a flav, energy fraction \a z and transverse
   * size \a r.
   */
  virtual Energy psi(int pol, int h, int hbar, int flav,
		     InvEnergy r, double z) = 0;

  /**
   * Return the square of the wavefunction summed over flavours and
   * helicities.
   */
  virtual Energy2 psi2(InvEnergy r, double z);
  //@}

public:

  /** @name Simple access functions. */
  //@{
  /**
   * Get the maxumim number of flavours considered. If negative only
   * the corresponding flavour will be considered.
   */
  inline int maxFlav() const;
  //@}

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
   * Called exactly once for each class by the class description
   * system before the main function starts or when this class is
   * dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();
  //@}

private:

  /**
   * The maxumim number of flavours considered. If negative only the
   * corresponding flavour will be considered.
   */
  int theMaxFlav;

protected:

  /**
   * The quark masses to be used (zero'th component is always ignored).
   */
  vector<Energy> qmass;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VectorMesonBase & operator=(const VectorMesonBase &);

};

}

#include "VectorMesonBase.icc"

#endif /* DIPSY_VectorMesonBase_H */
