// -*- C++ -*-
#ifndef PYTHIA7_StringRegion_H
#define PYTHIA7_StringRegion_H
// This is the declaration of the StringRegion class.

#include "FragConfig.h"
#include "String.h"
#include "Oriented.h"


namespace Pythia7 {

/**
 * The StringRegion describes a region of a String in the momentum
 * (parameter) space representation. During the fragmentation
 * procedure creation and destruction of <code>StringRegion</code>s is
 * the reponsibility of the String class.
 *
 * <b>Internal conventions</b>: They follow the <a
 * href="http://www.thep.lu.se/staff/torbjorn/Pythia.html">Pythia6</a>
 * conventions. Each StringRegion is ordered in the string
 * representation given its \f$(j,k)\f$ indices.  The \f$(j,k)\f$
 * <code>StringRegion</code> describes the region spanned by the
 * \f$(p_{+j},p_{-k})\f$ light-cone momenta.
 *
 * The forward, backward interfaces return the value of a
 * side-dependent variable given the side chosen to performe the
 * fragmentation (see String for conventions).
 * 
 */
class StringRegion {   

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline StringRegion(); 

  /**
   * Create the string region (\a j, \a k) of the string \a str,
   * spanned by the \f$p_{+j}\f$ and \f$p_{-k})\f$ light-cone momenta.
   */
  inline StringRegion(int j, int k, StringPtr str);

  /**
   * Destructor.
   */
  inline ~StringRegion();
  //@}

  /**
   * Return the four-vectors transverse to the string directions for 
   * this StringRegion
   */
  inline const LorentzVector<double> & ex() const;

  /**
   * Return the four-vectors transverse to the string directions for 
   * this StringRegion
   */
  inline const LorentzVector<double> & ey() const;

  /**
   * Return the Oriented forward light-cone momentum.
   */
  inline const LorentzMomentum & Pfwd() const;

  /**
   * Return the Oriented backward light-cone momentum.
   */
  inline const LorentzMomentum & Pbwd() const;

  /**
   * Return the energy squared available in this StringRegion.
   */
  inline Energy2 remW2() const;

  /**
   *Return the total energy squared of this StringRegion.
   */
  inline Energy2 W2() const;

  /**
   * Return the fraction of momentum remaning along the forward axis
   * of this StringRegion
   */
  inline double Xremf() const; 

  /**
   * Return the fraction of momentum remaning along the backward axis
   * of this StringRegion
   */
  inline double Xremb() const; 

  /**
   * Return the forward index.
   */
  inline int Ifwd() const;

  /**
   * Return the backward index.
   */
  inline int Ibwd() const;

  /**
   * Return the \f$j\f$ index of this StringRegion.
   */
  inline int j() const;

  /**
   * Return the \f$k\f$ index of this StringRegion.
   */
  inline int k() const;

  /**
   * Return true if this StringRegion is a primary (\f$j=k\f$) region.
   */
  inline bool aPrimaryStringRegion() const; 

  
  /**
   * Return the pointer to the String object.
   */
  inline StringPtr getTheStringPtr() const;

  /**
   * Update Function. (<b>Warning!</b> Temporary, not good!)
   */
  inline void setXrem(double , double ) const;
   
  /**
   * Print function used for debugging.
   */
  void echo() const;

private:

  /**
   * Initialize this region
   */
  void init(); 

  /**
   * Return the positive light-cone momentum of the StringRegion.
   */
  inline const LorentzMomentum& Pplus() const;

  /**
   * Return the negative light-cone momentum of the StringRegion.
   */
  inline const LorentzMomentum& Pminus() const;

private: 

  /**
   * The j index
   */
  int jj; 

  /**
   * The k index
   */
  int kk; 

  /**
   * The String
   */
  StringPtr theString;     

  /**
   *The transvers-Vector in x.
   */
  LorentzVector<double>  Ex;

  /**
   *The transvers-Vector in y.
   */
  LorentzVector<double>  Ey;

  /**
   *The StringRegion total energy squared.
   */
  Energy2 w2;                        

};

}

#include "StringRegion.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "StringRegion.tcc"
#endif

#endif /* PYTHIA7_StringRegion_H */
