// -*- C++ -*-
#ifndef ARIADNE5_ISGtoQEmission_H
#define ARIADNE5_ISGtoQEmission_H
//
// This is the declaration of the ISGtoQEmission class.
//

#include "ISQEmission.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ISGtoQEmission class contains all information about a generated
 * and performed initial-state gluon splitting into a q-qbar pair.
 */
class ISGtoQEmission: public ISQEmission {

public:

  /** @name Standard constructors and destructors. */
  //@{
   /**
   * The only relevant constructor.
   */
  ISGtoQEmission(const EmitterBase & inmodel, const DipoleBase & indipole,
		 tRemParPtr rem, tcPDPtr qin, Energy mqin,
		 const Lorentz5Momentum & ph)
    : ISQEmission(inmodel, indipole, rem, ph) {
    q = qin;
    mq = mqin;
    pold.first = rem->momentum();
    pold.second = ph;
  }

  /**
   * The deault constructor should not normally be used.
   */
  ISGtoQEmission() {}

  /**
   * The destructor.
   */
  virtual ~ISGtoQEmission();
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual ClonePtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual ClonePtr fullclone() const;
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
  ISGtoQEmission & operator=(const ISGtoQEmission &);

};

}

#endif /* ARIADNE5_ISGtoQEmission_H */
