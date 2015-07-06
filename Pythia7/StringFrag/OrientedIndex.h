// -*- C++ -*-
#ifndef PYTHIA7_OrientedIndex_H
#define PYTHIA7_OrientedIndex_H
// This is the declaration of the OrientedIndex class.

#include "Pythia7/Config/Pythia7.h"
// #include "OrientedIndex.fh"
// #include "OrientedIndex.xh"
#include "Oriented.h"



namespace Pythia7 {


/**
 * During the <i>stepping</i> procedure of the LundFragHandler one may
 * have to pass through several <code>StringRegion</code>s to step
 * from the initial StringRegion to the current one.  OrientedIndex is
 * an helper class used to iterate the indices \f$(j,k)\f$ in a
 * <i>side-independent</i> manner in order to allow for one generic
 * implementation of the stepping algorithmic.
 *
 * The OrientedIndex provides strait forward implementation of the
 * different operators needed for iteration, like <code>++</code>,
 * <code>--</code>, <code>>=</code>, <code>=< </code> etc. that use
 * the Oriented to access to the current side in the fragmentation.
 *
 * @see LundFragHandler
 * @see Oriented
 *
 */
class OrientedIndex: public Oriented{

public:

  /**
   * Standard constructor taking the start index as argument.
   */
  inline OrientedIndex(int ii=0): Idx(ii){}

  /**
   * Copu-constructor.
   */
  inline OrientedIndex(const OrientedIndex& OI): Idx(OI.Idx){}

  /**
   * Assignment operator.
   */
  inline void operator=(const OrientedIndex& OI){ 
    Idx = OI.Idx;
  }
  
  /**
   * Assignment operator for giving an index.
   */
  inline OrientedIndex& operator=(int ii) {
    Idx = ii;
    return *this;
  }

  /**
   * Increase or decrease index depending on Oriented.
   */
  inline int operator++() {
    return (Oriented::Dir()==Oriented::right)? ++Idx: --Idx;
  }

  /**
   * Increase or decrease index depending on Oriented.
   */
  inline int operator++(int) {
    return (Oriented::Dir()==Oriented::right)? Idx++: Idx--;
  }
  
  /**
   * Test for equality.
   */
  inline bool operator==(int ii) { 
    return (Idx ==ii);
  }

  /**
   * Test for inequality.
   */
  inline bool operator!=(int ii) const {
    return (Idx !=ii);
  }

  /**
   * Test for ordering.
   */
  inline bool operator < (int ii) const {
    return (Oriented::Dir()==Oriented::right)? 
      (Idx < ii) : (ii < Idx);
  }

  /**
   * Test for ordering.
   */
  inline bool operator <= (int ii) const {
    return (Oriented::Dir()==Oriented::right)?
      (Idx <= ii) : (ii <= Idx);
 }

  /**
   * Test for ordering.
   */
  inline bool operator < (const OrientedIndex& OI) const {
    return (Oriented::Dir()==Oriented::right)?
      (Idx < OI.idx()) : (OI.idx() < Idx);
  }

  /**
   * Test for ordering.
   */
  inline bool operator <= (const OrientedIndex& OI) const {
    return (Oriented::Dir()==Oriented::right)?
      (Idx <= OI.idx()) : (OI.idx() <= Idx);
  }

  /**
   * Negate the index.
   */
  inline OrientedIndex operator-() {
    return OrientedIndex(-Idx);
  }

  /**
   * Return the index.
   */
  inline int idx()const {
    return Idx;
  }

  /**
   * Return the index.
   */
  inline operator int() const{
    return idx();
  }

private :

  /**
   * The index.
   */
  int Idx;

};

}

//#include "OrientedIndex.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "OrientedIndex.tcc"
#endif

#endif /* PYTHIA7_OrientedIndex_H */
