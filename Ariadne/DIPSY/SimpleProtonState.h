// -*- C++ -*-
#ifndef DIPSY_SimpleProtonState_H
#define DIPSY_SimpleProtonState_H
//
// This is the declaration of the SimpleProtonState class.
//

#include "DipoleState.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the SimpleProtonState class.
 */
class SimpleProtonState: public DipoleState {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SimpleProtonState();

  /**
   * The standard constructor.
   */
  SimpleProtonState(const DipoleEventHandler & eh, Energy plus, Energy minus,
		    double angleWidth, double rapWidth, int connected,
		    WFInfoPtr wfi, double weight);

  /**
   * The copy constructor.
   */
  inline SimpleProtonState(const SimpleProtonState &);

  /**
   * The destructor.
   */
  virtual ~SimpleProtonState();
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
  SimpleProtonState & operator=(const SimpleProtonState &);

};

}

#include "SimpleProtonState.icc"

#endif /* DIPSY_SimpleProtonState_H */
