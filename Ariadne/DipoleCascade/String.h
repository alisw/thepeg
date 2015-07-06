// -*- C++ -*-
#ifndef ARIADNE_String_H
#define ARIADNE_String_H
//
// This is the declaration of the String class.
//

#include "CascadeBase.h"
#include "String.fh"
#include "Parton.fh"
#include "EMDipole.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * The String class respresentsa set of colour-connected Parton
 * objects. These are represented by the first and last Parton in the
 * string, where the first has anti-colour and the last has colour.
 */
class String: public CascadeBase {

public:

  /**
   * The DipoleState is a friend.
   */
  friend class DipoleState;

  /**
   * A pair of partons.
   */
  typedef pair<tParPtr,tParPtr> tParPair;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline String();

  /**
   * The copy constructor.
   */
  inline String(const String &);

  /**
   * The destructor.
   */
  inline virtual ~String();
  //@}

public:

  /**
   * Split the gluon, \a g, in this string into \a q and \a qbar. In
   * case this is not a closed string split the it into two. Also
   * connect the dipoles and partons correctly. The function also
   * removes the EM dipole connecting the endpoints of the string if
   * such a dipole exits. If \a remg is true the gluon will be from
   * the dipole state.
   */
  void split(tParPtr g, tParPtr q, tParPtr qbar, bool remg = false);

  /**
   * Move the endpoint of the string. The first argument is the old
   * endpoint and the second the new endpoint.
   */
  void moveEndpoint(tParPtr, tParPtr);

  /**
   * Join the current string to the string with the endpoint q and
   * replace the endpoints with the parton g. This function does not
   * remove the previous endpoints from the dipole state.
   */
  void join(tParPtr q, tParPtr g);

  /**
   * Return the EM dipole connected to the first and last parton
   * of the string. If none exists return NULL.
   */
  inline tEMDipPtr EMDip();

  /**
   * Set the EM dipole connected to the first and last parton
   * of the string.
   */
  inline void EMDip(tEMDipPtr x);

protected:

  /** @name Access the endpoints. */
  //@{
  /**
   * The first and last parton in this string. The first parton is
   * always anti-coloured, while the second is always coloured.
   */
  inline tParPair endpoints() const;

  /**
   * Set the \a first and \a last parton in this string. The first parton is
   * always anti-coloured, while the second is always coloured.
   */
  inline void endpoints(tParPtr first, tParPtr last);
  //@}
  
  /** @name Functions relating to the DipoleState to which this belongs. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual ClonePtr clone() const;

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  inline virtual void fillReferences(CloneSet &) const;

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
   * The first and last parton in this string. The first parton is
   * always anti-coloured, while the second is always coloured.
   */
  tParPair theEndpoints;

  /**
   * A pointer to the EM Dipole connected to the first and
   * last parton in the string if such a dipole exits.
   */
  tEMDipPtr theEMDipole;

public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<String> initString;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  String & operator=(const String &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of String. */
template <>
struct BaseClassTrait<Ariadne::String,1> {
  /** Typedef of the first base class of String. */
  typedef Ariadne::CascadeBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the String class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::String>
  : public ClassTraitsBase<Ariadne::String> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::String"; }
  /** Return the name of the shared library be loaded to get
   *  access to the String class and every other class it uses
   *  (except the base class). */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "String.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "String.tcc"
#endif

#endif /* ARIADNE_String_H */
