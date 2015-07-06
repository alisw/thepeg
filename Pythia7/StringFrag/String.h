// -*- C++ -*-
#ifndef PYTHIA7_String_H
#define PYTHIA7_String_H
// This is the declaration of the String class.

#include "FragConfig.h"
// #include "String.fh"
// #include "String.xh"
#include "ThePEG/EventRecord/Particle.h" 
#include "Oriented.h" 


namespace Pythia7 {

/** Two-dimensional index for <code>StringRegion</code>s. */
typedef pair<int,int> StringRegionIndex;
/** Map of <code>StringRegion</code>s keyed with pairs of integers. */
ThePEG_DECLARE_MAP(StringRegionIndex, StringRegionPtr, StringRegionMap);
  // typedef map<StringRegionIndex, StringRegionPtr> StringRegionMap;
/** Iterator for StringRegionMap. */
typedef StringRegionMap::iterator StrRgMapIt;


/**
 * The String class describes the string of the Lund string
 * fragmentation model.  implemented in the momentum (parameter) space
 * representation.
 *
 * The String describes both <i>open</i> of <i>closed</i> strings.
 * Given one end point particle of the string, String(cPPtr) uses the
 * colour-connection to build up the full string representation.
 * 
 * The string-regions are described by StringRegion objects.  They are
 * created dynamically when needed in the fragmentation process.  by
 * the getStringRegion() method.  When created for the first time, a
 * StringRegion is stored in a StringRegionMap according to the
 * indices \f$(j,k)\f$ of the region light-cone momentum pair
 * \f$(p_{+j} , p_{-k})\f$ available from the p_plus() and p_minus()
 * functions.
 *
 * Conventions for String construction:<br>
 * The PDPtr sent to to the string constructor defines the
 * <i>rightmost</i> endpoint of the string in the momentum space
 * representation. Then <code>Plus</code> denotes the direction
 * towards this <i>rightmost</i> end-point while <code>Minus</code>
 * denotes the opposite direction.<br>
 * String provides interfaces to access variables that depend on the
 * orientation chosen to perform a step in the fragmentation
 * procedure.  If the step is taken from the right of the plane
 * representation then a <code>fwd</code> (<code>bwd</code>) interface
 * to such a variable returns its <code>Plus</code>
 * (<code>Minus</code>) component, while if the step is taken from the
 * left then a <code>fwd</code> (<code>bwd</code>) interface returns
 * its <code>Minus</code> (<code>Plus</code>) component.<br>
 * 
 * <b> Warning :</b> No consistency checks are made by the String.  It
 * is the responsibility of the client using the String class to send
 * the correct end point to the string constructor.
 *
 * @see StringRegion
 * @see Oriented
 *
 */
class String {

  /** Vector of remaining momentum fractions. */
  typedef vector<double> xRemVec;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline String();

  /**
   * Copy-constructor.
   */
  inline String(const String &);

  /**
   * Standard constructor. Given the vector of particles of the string
   * computes the light-cone-like momenta \f$(p_+,p_-)\f$ of all
   * string pieces following the chain of the first end-point
   * anti-colour neighbours. Also provide a cutoff in invariant mass
   * squared below which two neighboring partons are clustered into
   * one.
   */
  String(const PVector &);

  /**
   * Destructor.
   */
  ~String();
  //@}

  /** @name Access momenta and momentum fractions. */
  //@{
  /**
   * Returns the total remaining momentum of this String, currently
   * available during fragmentation.
   */
  inline const LorentzMomentum&  PtotRem() const;

  /**
   * Returns the energy squared (Wrem2) of this String, currently
   * available during fragmentation.
   */
  inline Energy2 Wrem2() const;

  /**
   * Given a region index, \a fidx, return the forward longitudinal
   * light-cone momentum.
   */
  inline const LorentzMomentum& Pfwd(int fidx) const;
  /**
   * Given a region index, \a bidx, return the backward longitudinal
   * light-cone momentum
   */
  inline const LorentzMomentum& Pbwd(int bidx) const;

  /**
   * Given a region index \a fidx, return the fraction of the
   * longitudinal momentum available for fragmentation along the
   * forward axis of this region.
   */
  inline double Xremf(int fidx) const;              

  /**
   * Given a region index \a bidx return the fraction of the
   * longitudinal momentum available for fragmentation along the
   * backward axis of this region.
   */
  inline double Xremb(int bidx) const;

  /**
   * Given the string region index \a j, return the light-cone 
   * momenta \f$p_{+j}\f$.
   */
  inline const LorentzMomentum & Pplus(int j) const;

  /**
   * Given the string region index \a k, return the light-cone 
   * momenta \f$p_{-k}\f$.
   */
  inline const LorentzMomentum & Pminus(int k) const;

  /**
   * Given the indicex \a j, return the remaining fraction of the
   * light-cone momenta \f$p_{+j}\f$ in the corresponding region.
   */
  inline double XplusRem(int j) const;

  /**
   * Given the indicex \a k, return the remaining fraction of the
   * light-cone momenta \f$p_{-k}\f$ in the corresponding region.
   */
  inline double XminusRem(int k) const;
  //@}

  /** @name Access <code>StringRegion</code>s. */
  //@{
  /**
   * Return the number of string pieces in this String
   */
  inline int nPrimaryStringRegion() const;

  /**
   * Given the indices \a j and \a k, return a pointer to the
   * corresponding StringRegion.
   */
  cStringRegionPtr getStringRegionPtr(int j, int k);

  /**
   * Given a pointer to a StringRegion (\a firstSR) return a pointer
   * to its next nearest StringRegion in the <i>up</i> direction.
   */
  cStringRegionPtr nextUp(cStringRegionPtr firstSR);

  /**
   * Given a pointer to a StringRegion (\a firstSR) return a pointer
   * to its next nearest StringRegion in the <i>down</i> direction.
   */
  cStringRegionPtr nextDown(cStringRegionPtr firstSR);

  /**
   * Return pointer to the first (rightmost) StringRegion in the
   * String plane representation.
   */
  inline cStringRegionPtr firstSR();  

  /**
   * Return pointer to the last (leftmost) StringRegion in the String
   * plane representation.
   */
  inline cStringRegionPtr lastSR();
  //@}

  /** @name Functions for updating the string (regions). */
  //@{
  /**
   * Compute thr string total momentum and initialize all Xrem(+,-)[i] to 1
   */
  inline void reset();

  /**
   * Compute the remaining string total momentum, each time a new hadron 
   * of momentum HadronP is produced in the fragmentation process.
   */
  inline void updatePtotrem(const LorentzMomentum& HadronP);

  /**
   * Given a hadron produced with fraction \a xfnH of the forward
   * momentum, compute the remaining fraction of the corresponding
   * string regions with forward index \a fidx.
   */
  inline void updateXremf(int fidx, double xfnH);

  /**
   * Given a hadron produced with fraction \a xbnH of the backward
   * momentum, compute the remaining fraction of the corresponding
   * string regions with backward index \a bidx.
   */
  inline void updateXremb(int bidx, double xbnH);

  /**
   * Set the remaining forward fraction of regions with index \a jj to
   * \a newX.  (Only used for closed string)
   */
  inline void setXplusRem(int jj, double newX);

  /**
   * Set the remaining backward fraction of regions with index \a kk
   * to \a newX.  (Only used for closed string)
   */
  inline void setXminusRem(int kk, double newX);
  
  //@}

  /**
   * Print function to be used for debugging.
   */
  void echo() const; 


private: 

  /**
   * Given the pointer to the first end-point computes the total  
   * momentum  of this String.
   */
  //  inline void sumMomentun(cPPtr firstPPtr);

  /**
   *Initialize all Xrem(i) fractions to 1
   */
  inline void InitXrem();

  /**
   * Delete all StringRegion objects and pointers of the StringRegionMap.
   */
  void clearStringRegionMap();

private:

  /**
   * The Total momentum of the String.
   */
  LorentzMomentum Ptotrem;

  /**
   * Vector of the light-cone momenta \f$p_{+j}\f$.     
   */
  MomentumVector pplus;

  /**
   * Vector of the light-cone momenta \f$p_{-k}\f$.     
   */
  MomentumVector pminus;      

  /**
   * Vector of the remaining positive light-cone fractions. 
   */
  xRemVec xPlusRem;

  /**
   * Vector of the remaining positive light-cone fractions. 
   */
  xRemVec xMinusRem;   

  /**
   * The map of <code>StringRegion</code>s.
   */
  StringRegionMap theStringRegionMap;

  /**
   * Number of primary string regions in this String.
   */
  int ns;

  /**
   * The initial total momentum.
   */
  LorentzMomentum Ptot;

  /**
   * The initial, possibly pre-clustered particles.
   */
  PVector partons;

};

}

#include "String.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "String.tcc"
#endif

#endif /* PYTHIA7_String_H */
