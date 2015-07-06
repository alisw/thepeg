// -*- C++ -*-
#ifndef PYTHIA7_EndPoint_H
#define PYTHIA7_EndPoint_H
// This is the declaration of the EndPoint class.

#include "FragConfig.h"
#include "ThePEG/PDT/ParticleData.h"
#include "StringRegion.h"
// #include "EndPoint.fh"
// #include "EndPoint.xh"

namespace Pythia7 {

/**
 * The EndPoint is the main helper class used by the LundFragHandler
 * to handle the steps in the fragmentation procedure. The endpoint
 * describes the quark (diquark) of a quark-antiquark
 * (diquark-antidiquark) pair created at a breakpoint during the
 * string fragmentation. At each step a new Hadron is created joining
 * the newly created EndPoint with the <i>old</i> EndPoint left over
 * at the previous step in the iterative procedure.
 * 
 * The information about the particle type are accessed through a
 * pointer to the ParticleData object.
 * 
 * Conventions of the stepping methods :<br>
 * stepUp() and stepDown() are the Oriented methods to move an
 * EndPoint through the grid of the StringRegions.  Internally the
 * EndPoint's StringRegion is defined in the grid by its
 * <tt>(j,k)</tt> indices.  Assume an EndPoint <tt>ep</tt> standing in
 * the <tt>(j<sub>i</sub>, k<sub>i</sub>)</tt>:<br>
 * then if a step is taken from the <i>right side</i> of the string :<br>
 * <code>ep.stepUp()</code> will move it to the region
 * <tt>(j<sub>i</sub>, k<sub>i+1</sub>)</tt> while
 * <code>ep.stepDown()</code> to the region <tt>(j<sub>i+1</sub>,
 * k<sub>i</sub>)</tt>.<br>
 * If a step is taken from the <i>left side</i>:<br>
 * <code>ep.stepUp()</code> will now move it to the region
 * <tt>(j<sub>i-1</sub>, k<sub>i</sub>)</tt> and
 * <code>ep.stepDown()</code> to the region <tt>(j<sub>i</sub>,
 * k<sub>i-1</sub>)</tt> <br>
 *
 * Apart of the pointer to the StringRegion, the step methods do not
 * change any other members of the EndPoint.
 *
 *
 * @see String
 * @see StringRegion
 * @see LundFragHandler
 * @see Oriented
 *
 */
class EndPoint{

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline EndPoint(); 

  /**
   * The copy constructor.
   */
  inline EndPoint(const EndPoint&);  

  /**
   * Create an EndPoint given the type () of the particle standing 
   * at this endpoint and the string region of the break point 
   */
  inline EndPoint(cPDPtr inPDPtr, cStringRegionPtr inSRPtr);     

  /**
   * The destructor.
   */
  inline ~EndPoint();

  /**
   * Reset this EndPoint to the default values.
   */
  inline void Init();

  /**
   * Reset this EndPoint's ParticleData and StringRegion.
   */
  inline void Init(cPDPtr inPDPtr, cStringRegionPtr inSRPtr);     

  /**
   * Return the charge conjugate of this EndPoint.
   */
  EndPoint CC() const; 

  /**
   * Stepping Method to move the EndPoint particle through the 
   * different string regions of the String representation
   */
  inline void stepUp();

  /**
   * Stepping Methods to move the EndPoint particle through the 
   * different string regions of the String representation
   */
  inline void stepDown();

  /**
   * Update this EndPoint given the current EndPoint (currEP). [At the
   * end of a step the current EndPoint becomes the last EndPoint for
   * the next step.]
   */
  void UpdatedFrom(const EndPoint& currEP);
 
  /**
   * Get the pointer to the StringRegion where this EndPoint is standing
   * during the fragmentation process
   */
  inline cStringRegionPtr SR() const; 
 
  /**
   * Set the pointer to the StringRegion where this EndPoint is standing
   * during the fragmentation process
   */
  inline void SR(cStringRegionPtr newSR);

  /**
   * Get the forward light-cone vector of the current StringRegion of
   * this EndPoint.
   */
  inline const LorentzMomentum & Pfwd() const;

  /**
   * Get the backward light-cone vector of the current StringRegion of
   * this EndPoint.
   */
  inline const LorentzMomentum & Pbwd() const;

  /**
   * Get the remaining forward momentum fraction of the current
   * StringRegion of this EndPoint
   */
  inline double Xremf() const;

  /**
   * Get the remaining backward momentum fraction of the current
   * StringRegion of this EndPoint
   */
  inline double Xremb() const;

  /**
   * Set the ParticleData pointer at this EndPoint
   */
  inline void PData(cPDPtr pd);

  /**
   * Get the ParticleData pointer at this EndPoint
   */
  inline cPDPtr PData()const;

  /**
   * Return the <em>constituent</em> mass of the Particle at this
   * EndPoint
   */
  inline Energy mass() const;
  
  /**
   * Set the transverse momentum components of the Particle at this
   * EndPoint
   */
  inline void setPt(const TransverseMomentum & newPT);

  /**
   * Set the transverse momentum components of the Particle at this
   * EndPoint
   */
  inline void setPt(Energy px, Energy py);

  /**
   * Get the transverse momentum components of the Particle at this
   * EndPoint
   */
  inline const TransverseMomentum &  pTcomp() const;

  /**
   * Get the transverse momentum components of the Particle at this
   * EndPoint
   */
  inline Energy Px() const;

  /**
   * Get the transverse momentum components of the Particle at this
   * EndPoint
   */
  inline Energy Py() const;

  /**
   * Return the transverse 4-vector of the endpoint particle in 
   * the current string region 
   */
  inline LorentzMomentum pT() const;

  /**
   * Set the \f$\Gamma\f$ variable of the endpoint
   * vertex. \f$\Gamma\f$ represents the invariant-time squared from
   * the vertex of the old EndPoint to this one
   */
  inline void Gamma(Energy2 );

  /**
   * Get the \f$\Gamma\f$ variable of the endpoint
   * vertex. \f$\Gamma\f$ represents the invariant-time squared from
   * the vertex of the old EndPoint to this one
   */
  inline Energy2 Gamma() const;

  /**
   * EndPoint assigment operator
   */
  inline EndPoint& operator = (const EndPoint& ep);


  /**
   * Print info to cout for debugging purposes.
   */
  void echo() const;

private:

  /**
   * Get the pointer to the string being fragmented
   */
  inline StringPtr theString();

private:

  /**
   * The particle type at this endpoint
   */
  cPDPtr theParticle;

  /**
   * The current string region 
   */
  cStringRegionPtr theStrRg;

  /**
   * The transverse momentum components of the particle  
   */
  TransverseMomentum thePTcomp;

  /**
   *The invariant-time square
   */
  Energy2 theGamma;

};


}

#include "EndPoint.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "EndPoint.tcc"
#endif

#endif /* PYTHIA7_EndPoint_H */









