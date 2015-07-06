// -*- C++ -*-
#ifndef Ariadne5_ColourChargeRegions_H
#define Ariadne5_ColourChargeRegions_H
//
// This is the declaration of the ColourChargeRegions class.
//

#include "Ariadne/Cascade/ReweightBase.h"
#include "Ariadne/Cascade/QCDDipole.fh"
#include "Ariadne/Cascade/Emission.fh"
namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ColourChargeRegions class can be used to reweight gluon
 * emissions so that a more correct colour charge is assigned to an
 * emission probability, based not only on the charge of the emitting
 * partons in the dipole, but also on where in rapidity the gluon is
 * emitted. This uses a simplified version of the model by Eden and
 * Gustafson [hep-ph/9805228].
 *
 * @see \ref ColourChargeRegionsInterfaces "The interfaces"
 * defined for ColourChargeRegions.
 */
class ColourChargeRegions: public Ariadne5::ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ColourChargeRegions();

  /**
   * The destructor.
   */
  virtual ~ColourChargeRegions();
  //@}

public:

  /** @name Virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * If the emission is from a q-qbar dipole, preweight with Nc/2/CF.
   */
  virtual double preweight(const Emission & emission) const;
 
  /**
   * Reweight with 2CF/Nc if the repisity of the gluon corresponds to
   * it being emitted off a quark line.
   */
  virtual double reweight(const Emission & emission) const;
  //@}

protected:

  /**
   * Set up the charge regions necessary to calculate the colour charge.
   */
  void setup(const QCDDipole & d) const;

  /**
   * Return the colour charge for the given rapidity.
   */
  double charge(double y) const;

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
   * The unique id of last emission for which the colour regions were
   * set up.
   */
  mutable cEmPtr lastEmission;

  /**
   * The map tracing the history of the coloured Parton in the Dipole.
   */
  mutable map<double,double> colregions;

  /**
   * The map tracing the history of the anti-coloured Parton in the Dipole.
   */
  mutable map<double,double> acoregions;

public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ColourChargeRegions & operator=(const ColourChargeRegions &);

};

}

#endif /* Ariadne5_ColourChargeRegions_H */
