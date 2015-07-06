// -*- C++ -*-
#ifndef ARIADNE_CascadeBase_H
#define ARIADNE_CascadeBase_H
//
// This is the declaration of the CascadeBase class.
//

#include "Ariadne/Config/CloneBase.h"
#include "CascadeBase.fh"
#include "DipoleState.fh"
#include "Ariadne/DipoleCascade/CascadeHandler.h"

namespace Ariadne {

using namespace ThePEG;

/**
 * CascadeBase is the base class of all Partons, Emitters and
 * DipoleState classes in the Ariadne dipole cascade classes. It keeps
 * track of the DioleState to which an object belongs, and has also a
 * pointer to the AriadneCascade object administering the cascade.
 */
class CascadeBase: public CloneBase {

public:

  /**
   * The DipoleState is a friend.
   */
  friend class DipoleState;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor has an optional pointer to the
   * CascadeHandler in charge as argument.
   */
  inline CascadeBase(tHandlerPtr hdl = tHandlerPtr());

  /**
   * The copy constructor.
   */
  inline CascadeBase(const CascadeBase &);

  /**
   * The destructor.
   */
  inline virtual ~CascadeBase();
  //@}

protected:

  /** @name Functions relating to the DipoleState and CascadeHandler
   *  to which this belongs. */
  //@{
  /**
   * Get the Ariadne::CascadeHandler in charge of the current generation.
   */
  inline tHandlerPtr handler() const;

  /**
   * Set the Ariadne::CascadeHandler in charge of the current generation.
   */
  inline void handler(tHandlerPtr);

  /**
   * Get the DipoleState to which this Dipole belongs.
   */
  inline tDipoleStatePtr state() const;

  /**
   * Set the DipoleState to which this Dipole belongs.
   */
  inline void state(tDipoleStatePtr);

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
  //@}

public:

  /** @name Functions determining if the object has been changed since
      the last generation. */
  //@{
  /**
   * If true, this object has been modified since the last round of
   * generating emissions.
   */
  inline bool touched() const;

  /**
   * Signal that this object has been modified.
   */
  inline void touch();

  /**
   * Signal that all possible emissions involving this object has been
   * generated.
   */
  inline void untouch();
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
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The Ariadne::CascadeHandler in charge of the current generation.
   */
  tHandlerPtr theHandler;

  /**
   * The DipoleState to which this Dipole belongs.
   */
  tDipoleStatePtr theState;

  /**
   * If true, this object has been modified since the last round of
   * generating emissions.
   */
  bool isTouched;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<CascadeBase> initCascadeBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CascadeBase & operator=(const CascadeBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CascadeBase. */
template <>
struct BaseClassTrait<Ariadne::CascadeBase,1> {
  /** Typedef of the first base class of CascadeBase. */
  typedef Ariadne::CloneBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CascadeBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::CascadeBase>
  : public ClassTraitsBase<Ariadne::CascadeBase> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::CascadeBase"; }
  /** Return the name of the shared library be loaded to get
   *  access to the CascadeBase class and every other class it uses
   *  (except the base class). */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "CascadeBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CascadeBase.tcc"
#endif

#endif /* ARIADNE_CascadeBase_H */
