// -*- C++ -*-
//
// XSecStat.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_XSecStat_H
#define THEPEG_XSecStat_H
//
// This is the declaration of the XSecStat class.
//

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * XSecStat is a concrete helper class used to collect statistics
 * about the cross section for a specific process or group of
 * processes. It contains an overestimated cross section and
 * information about the number of times the process has been used to
 * generate an event and how many times this event has been accepted.
 *
 * An object of this class must initially be given an overestimated
 * cross section in the constructor or with the maxXSec(CrossSection)
 * function. Each time the corresponding process is selected
 * (according to maxXSec()), the select(double) function should be
 * called giving the weight with which the event will be accepted as
 * argument. If the event is then accepted, the accept() function
 * should be called. If an event is later vetoed, the reject()
 * function should be called.
 * 
 */
class XSecStat {

public:

  /**
   * Enumerate the different weight classes
   */
  enum {
    plainWeights = 0,
    plainVetoedWeights,
    reweightedWeights,
    reweightedVetoedWeights
  };

  /** @name Standard constructors, destructor and assignment operator. */
  //@{
  /**
   * The default constructor.
   */
  XSecStat() 
    : theMaxXSec(ZERO), theAttempts(0), theAccepted(0), theVetoed(0),
      theSumWeights (),
      theSumWeights2(), theLastWeight(0.0) {}

  /**
   * Constructor taking the overestimated cross section, \a xsecmax,
   * as argument.
   */
  explicit XSecStat(CrossSection xsecmax) 
    : theMaxXSec(xsecmax), theAttempts(0), theAccepted(0), theVetoed(0),
      theSumWeights (),
      theSumWeights2(), theLastWeight(0.0) {}

  /**
   * The assignment operator.
   */
  XSecStat & operator=(const XSecStat & x) = default;

  /**
   * Add the contents of another XSecStat.
   */
  XSecStat & operator+=(const XSecStat & x) {
    theAttempts    += x.theAttempts;
    theAccepted    += x.theAccepted;
    theVetoed      += x.theVetoed;
    for( unsigned int ix = 0; ix < 4; ++ix ) {
      theSumWeights [ix] += x.theSumWeights [ix];
      theSumWeights2[ix] += x.theSumWeights2[ix];
    }
    theLastWeight = 0.0;
    return *this;
  }

  /**
   * Reset the statistics.
   */
  void reset() {
    theAttempts = theAccepted = theVetoed = 0;
    theSumWeights = theSumWeights2 = {};
    theLastWeight = 0.0;
  }

  //@}

public:

  /** @name Simple access functions */
  //@{

  /**
   * An event of the corresponding class has been accepted. The
   * select() method must have been called before.
   */
  void accept() { 
    theAccepted += 1;
  }

  /**
   * An event of the corresponding class has been attempted. It will
   * subsequently be accepted with the given \a weight.
   */
  void select(double weight) {
    theAttempts += 1;
    theSumWeights [reweightedWeights] +=     weight ;
    theSumWeights2[reweightedWeights] += sqr(weight);
    theSumWeights [plainWeights]      +=     weight ;
    theSumWeights2[plainWeights]      += sqr(weight);
    theLastWeight = weight;
  }

  /**
   * Reweight a selected and accepted event.
   */
  void reweight(double oldWeight, double newWeight) {
    theSumWeights [reweightedWeights] +=     newWeight  -     oldWeight ;
    theSumWeights2[reweightedWeights] += sqr(newWeight) - sqr(oldWeight);
  }

  /**
   * Reject the event which was last accepted with accept() or
   * selected with select(double). The \a weight should be set to the
   * value, \f$w\f$, used in the previous call to select(double),
   * except if the event has been accepted with the probability
   * \f$w\f$, in which case \a weight should be set to \f$sign(1,
   * w)\f$.
   */
  void reject(double weight = 1.0) {
    theSumWeights [reweightedVetoedWeights] +=            weight ;
    theSumWeights2[reweightedVetoedWeights] +=        sqr(weight);
    theSumWeights [plainVetoedWeights]      +=     theLastWeight ;
    theSumWeights2[plainVetoedWeights]      += sqr(theLastWeight);
    theVetoed += 1;
  }

  /**
   * The overestimated cross section.
   */
  CrossSection maxXSec() const { return theMaxXSec; }

  /**
   * The sum of the weights so far.
   */
  double sumWeights() const { 
    return theSumWeights[reweightedWeights] - theSumWeights[reweightedVetoedWeights];
  }

  /**
   * The sum of the squared weights so far.
   */
  double sumWeights2() const {
    return theSumWeights2[reweightedWeights] + theSumWeights2[reweightedVetoedWeights];
  }

  /**
   * The sum of the weights so far, excluding reweighting.
   */
  double sumWeightsNoReweight() const { 
    return theSumWeights[plainWeights] - theSumWeights[plainVetoedWeights];
  }

  /**
   * The sum of the squared weights so far, excluding reweighting.
   */
  double sumWeights2NoReweight() const { 
    return theSumWeights2[plainWeights] + theSumWeights2[plainVetoedWeights];
  }

  /**
   * The current estimate of the cross section for the corresponding
   * class of events. If no events have been generated, maxXSec() will
   * be returned.
   */
  CrossSection xSec(double att = 0) const {
    double n = (att == 0.0 ? attempts() : att);
    return n ? maxXSec()*sumWeights()/n : maxXSec();
  }

  /**
   * The current estimate of the error in the cross section for the
   * corresponding class of events. If no events have been generated,
   * maxXSec() will be returned.
   */
  CrossSection xSecErr(double att = 0) const {
    double n = (att == 0.0 ? attempts() : att);
    if ( n < 2 )
      return maxXSec();
    double sw = sumWeights(); double sw2 = sumWeights2();
    return
      maxXSec()*sqrt(abs(sw2/n-sqr(sw/n))/(n-1));
  }

  /**
   * The current estimate of the cross section for the corresponding
   * class of events, excluding reweighting. If no events have been
   * generated, maxXSec() will be returned.
   */
  CrossSection xSecNoReweight(double att = 0) const {
    double n = (att == 0.0 ? attempts() : att);
    return n ? maxXSec()*sumWeightsNoReweight()/n : maxXSec();
  }

  /**
   * The current estimate of the error in the cross section for the
   * corresponding class of events, excluding reweighting. If no
   * events have been generated, maxXSec() will be returned.
   */
  CrossSection xSecErrNoReweight(double att = 0) const {
    double n = (att == 0.0 ? attempts() : att);
    if ( n < 2 )
      return maxXSec();
    double sw = sumWeightsNoReweight(); 
    double sw2 = sumWeights2NoReweight();
    return
      maxXSec()*sqrt(abs(sw2/n-sqr(sw/n))/(n-1));
  }

  /**
   * Number of attempts so far.
   */
  double attempts() const { return theAttempts; }

  /**
   * Number of accepts so far.
   */
  double accepted() const { return theAccepted-theVetoed; }

  /**
   * Number of vetoes so far.
   */
  double vetoed() const { return theVetoed; }

  /**
   * Set the overestimated cross section.
   */
  void maxXSec(CrossSection x) { theMaxXSec = x; }
  //@}

public:

  /** @name I/O functions */
  //@{
  /**
   * Output to a persistent stream.
   */
  void output(PersistentOStream & os) const;

  /**
   * Input from a persistent stream.
   */
  void input(PersistentIStream & is);
  //@}

private:

  /**
   * The overestimated cross section.
   */
  CrossSection theMaxXSec;

  /**
   * Number of attempts so far.
   */
  double theAttempts;

  /**
   * Number of accepted events so far.
   */
  double theAccepted;

  /**
   * Number of events vetoed after being accepted
   */
  double theVetoed;

  /**
   * The sum of the weights so far.
   */
  array<double,4> theSumWeights;

  /**
   * The sum of the squared weights so far.
   */
  array<double,4> theSumWeights2;

  /**
   * The last selected weight, ignoring reweighting.
   */
  double theLastWeight;

};

/** Ouptut an XSecStat to a persistent stream. */
PersistentOStream & operator<<(PersistentOStream &, const XSecStat &);

/** Input an XSecStat from a persistent stream. */
PersistentIStream & operator>>(PersistentIStream &, XSecStat &);

/** Add the contents of two XSecStat objects. */
inline XSecStat operator+(const XSecStat & x1, const XSecStat & x2) {
  XSecStat x = x1;
  return x += x2;
}

}

#endif /* THEPEG_XSecStat_H */
