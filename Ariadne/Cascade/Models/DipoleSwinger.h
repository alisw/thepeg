// -*- C++ -*-
#ifndef ARIADNE5_DipoleSwinger_H
#define ARIADNE5_DipoleSwinger_H
//
// This is the declaration of the DipoleSwinger class.
//

#include "Ariadne/Cascade/EmitterBase.h"
#include "DipoleSwing.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/PartonTraits.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The DipoleSwinger class implements the final-state swing method for
 * colour reconections.
 *
 * @see \ref DipoleSwingerInterfaces "The interfaces"
 * defined for DipoleSwinger.
 */
class DipoleSwinger: public EmitterBase {

public:

  /**
   * Convenient typedef.
   */
  ThePEG_DECLARE_POINTERS(Ariadne5::DipoleSwing,DipSwPtr);

  /**
   * Another convenient typedef.
   */
  typedef map<tcQCDPtr,DipSwPtr> CacheMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleSwinger();

  /**
   * The destructor.
   */
  virtual ~DipoleSwinger();
  //@}

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if and only if this model can handle the given
   * Emitter.
   */
  virtual bool canHandle(const DipoleBase &) const;

  /**
   * If the given EmissionModel overlaps with this model for the given
   * Emitter, return true if this model should take precedence. Must
   * only be called for Emitters for which canHandle() is true.
   */
  virtual bool overrides(const EmitterBase &, DipoleBase &) const;

  /**
   * Check if objects related to the given \a dipole have been touched
   * in a way such that emissions must be regenerated.
   */
  virtual bool touched(const DipoleBase & dipole) const;

  /**
   * Generate the a phase space point for an emission corresponding to
   * this model. Must only be called for Emitters for which
   * canHandle() is true.
   */
  virtual EmPtr generate(const DipoleBase &, Energy rhomin, Energy rhomax) const;

  /**
   * Perform an emission previously generated for this Emitter. Must
   * only be called for Emitters for which canHandle() is true.
   * @return true if the emission was successful
   */
  virtual bool perform(const Emission &) const;

  /**
   * Reverse a previously performed emission. Sub-classes which has
   * signalled that they can revert an emission but fails to do so,
   * must throw a Exception::runerror.
   */
  virtual void revert(const Emission & emission) const;

  /**
   * Possibility to override the general non-perturbative cutoff for
   * this process. The default is AriadneHandler::pTCut() which is
   * typically used for QCD radiation.
   */
  virtual Energy rhoCut() const;

  /**
   * The maximum allowed value of the evolution variable, Used to
   * disallow perturbative swings. Inactive if below rhoCut().
   */
  inline Energy maxRho() const {
    return theMaxRho;
  }

  //@}

  /**
   * Swing the given dipoles.
   */
  static void swing(tQCDPtr d1, tQCDPtr d2);

  /**
   * Find strings with total mass below some cutoff and place them in
   * the tinyString container.
   */
  void findTinyStrings() const;

public:

  /**
   * Select the given Swing if it has the largerst scale.
   */
  void select(DipSwPtr sel) const {
    if ( sel && ( !lastSelected || sel->rho > lastSelected->rho ) ) lastSelected = sel;
  }


  /**
   * Internal helper class to keep track of distnaces between partons.
   */
  struct DeltaTau2 {

    /**
     * Constructor taked two partons and calculates the variables
     * needed to get the distance as a function of time.
     */
    DeltaTau2(tcParPtr p1, tcParPtr p2, int optin, Length Rmax): sizeOpt(optin) {
      LorentzDistance dx = p1->vertex() - p2->vertex();
      dx2 = min(dx.m2(), ZERO);
      LorentzVector<double> dp =
	p1->momentum()/p1->momentum().e() - p2->momentum()/p2->momentum().e();
      dp2 = dp.m2();
      dpdotdx = dp*dx;
      if ( sizeOpt >= 1 )
	dp2 = max((p1->momentum() + p2->momentum()).m2(), sqr(hbarc/Rmax))/GeV2;
      pp1 = p1;
      pp2 = p2;
    }

    /**
     * Return distance between partons after a given time \a dt.
     */
    Area dTau2(Time dt) const {
      if ( sizeOpt <= 0 )
	return - dx2 - dt*2.0*dpdotdx - sqr(dt)*dp2;
      else
	return sqr(dt)*dp2;
    }

    /**
     * Return the impact-parameter distance between partons after a
     * given time \a dt.
     */
    Length db() const {
      return sqrt(-dx2);
    }

    /**
     * Return the impact-parameter distance between partons after a
     * given time \a dt.
     */
    Area db2() const {
      return -dx2;
    }

    /**
     * Function version of dTau2
     */
    Area operator()(Time dt) const {
      return dTau2(dt);
    }

    /**
     * If the distance between the partons may become zero in the
     * given time interval, return the corresponding times. If there
     * are no zeros, time=0 will be returned.
     */
    pair<Length,Length> zeros(Time dt1, Time dt2) const {
      pair<Length,Length> ret;
      Area squarg = sqr(dpdotdx) - dx2*dp2;
      if ( squarg < ZERO ) return ret;
      Time t0 = - (dpdotdx + sqrt(squarg))/dp2;
      if ( dt1 <= t0 && t0 <= dt2 ) ret.first = t0;
      t0 = - (dpdotdx - sqrt(squarg))/dp2;
      if ( dt1 <= t0 && t0 <= dt2 ) ret.second = t0;
      
      return  ret;
    }

    /**
     * Return the minimum and maximum in the given interval.
     */
    pair<Area,Area> minmax(Time dt1, Time dt2) const {
      pair<Area,Area> ret(dTau2(dt1),dTau2(dt2));
      if ( ret.first > ret.second ) swap(ret.first, ret.second);
      Time dtext = -dpdotdx/dp2;
      if ( dtext >= dt2 || dt1 >= dtext ) return ret;
      Area dtau2ext = dTau2(dtext);
      if ( dtau2ext > ret.second ) ret.second = dtau2ext;
      if ( dtau2ext < ret.first ) ret.first = dtau2ext;
      return ret;
    }

    /**
     * Return the extreme point if in the given interval.
     */
    Time extreme(Time dt1, Time dt2) const {
      return min(max(dt1, -dpdotdx/dp2), dt2);
    }

    /**
     * The initial distance between the partons.
     */
    Area dx2;

    /**
     * The distance between the momentum vectors, scaled by their energies. 
     */
    double dp2;

    /**
     * Different options for calculating the dipole size.
     */
    int sizeOpt;

    /**
     * The scalar product between the difference between momentum
     * vectors (scaled with their energies) and the difference between
     * the initial positions.
     */
    Length dpdotdx;

    /* **** ATTENTION *** Only for debugging */
    InvEnergy2 chris(Time dt) const {
      InvEnergy2 time = sqr(dt/hbarc);
      Energy2 E2 = pp1->momentum().plus()*pp2->momentum().minus() +
	           pp1->momentum().minus()*pp2->momentum().plus();
      InvEnergy x1x = pp1->vertex().x()/hbarc -  pp2->vertex().x()/hbarc + time*pp1->momentum().x();
      InvEnergy x1y = pp1->vertex().y()/hbarc -  pp2->vertex().y()/hbarc + time*pp1->momentum().y();
      InvEnergy x2x = pp1->vertex().x()/hbarc -  pp2->vertex().x()/hbarc - time*pp2->momentum().x();
      InvEnergy x2y = pp1->vertex().y()/hbarc -  pp2->vertex().y()/hbarc - time*pp2->momentum().y();
      InvEnergy2 xxScalar = x1x*x2x + x1y*x2y;
      Energy2 ppScalar = pp1->momentum().x()*pp2->momentum().x() +
	                 pp1->momentum().y()*pp2->momentum().y();
      return sqr(time)*E2 + xxScalar - sqr(time)*ppScalar;
    }

    tcParPtr pp1, pp2;

  };

  /**
   * Internal helper class to keep track of ratio of distances between
   * original and reconected dipoles.
   */
  struct TauRatio {

    /**
     * The constructor takes the original dipole pair as argument.
     */
    TauRatio(const QCDDipole & d1, const QCDDipole & d2, Length Rmaxin, int optin);

    /**
     * Return the swing ratio at the given time.
     */
    double rat(Time dt) const;

    /**
     * Return the swing ratio at the given time.
     */
    double operator()(Time dt) const {
      return rat(dt);
    }

    /**
     * Estimate a maximum value of the swing ratio.
     */
    double maximum(Time dt1, Time dt2) const;

    void plot(Time dt1, Time dt2, int N = 20) const;

    /**
     * The objects for calculating distances between the original and
     * reconnected dipoles.
     */
    DeltaTau2 d1tau2, d2tau2, r1tau2, r2tau2;

    /**
     * Hadronic size to regularize large dipoles.
     */
    Length Rmax;

    /**
     * Different options for calculating the dipole size.
     */
    int sizeOpt;

    /* **** ATTENTION *** Only for debugging */
    void debug() const;
    double chris(Time dt) const {
      return d1tau2.chris(dt)*d2tau2.chris(dt)/(r1tau2.chris(dt)*r2tau2.chris(dt));
    }
    double srat;
    Energy2 s12, s34, s14, s32;

  };


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

protected:

  /**
   * The parameter giving the rate of swinging
   */
  double lambda;

  /**
   * The confinement scale used to prevent swinging to too large dipoles.
   */
  Length Rmax;

  /**
   * Should we use linear or logarithmic evolution in time?
   */
  bool linear;

  /**
   * Different options for calculating the dipole size.
   */
  int sizeOpt;

  /**
   * Cutoff in evolution the variable, if different from the standard QCD cutoff.
   */
  Energy theRhoCut;

  /**
   * The maximum allowed value of the evolution variable, Used to
   * disallow perturbative swings. Inactive if below rhoCut().
   */
  Energy theMaxRho;

  /**
   * Cutoff against too small string. No swing is allowed to produce
   * strings with invarint mass below this value.
   */
  Energy stringMassCut;

  /**
   * The tiniest string in the dipole system if below stringMassCut.
   */
  mutable set<tcDBPtr> tinyString;

  /**
   * Previously calculated cached swings.
   */
  mutable CacheMap cache;

  /**
   * A pointer to the DipoleState last treated.
   */
  mutable tDipoleStatePtr lastState;

  /**
   * The uniqueId of DipoleState last treated.
   */
  mutable unsigned long lastStateId;

  /**
   * A pointer to the last selected generated swing.
   */
  mutable DipSwPtr lastSelected;

  /**
   * A list of swing dipoles that have been rejected and should
   * therefore not be tried again until something has changed in the
   * dipole state.
   */
  mutable set< pair<tQCDPtr,tQCDPtr> > disallowed;

  /**
   * If something has been entered in the disallowed set of dipole
   * pairs, this was the maximum scale when they were entered,
   */
  mutable Energy disallowedrho;

private:

  /**
   * Helper function used by the interface.
   */
  string setRmax(string);

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleSwinger & operator=(const DipoleSwinger &);

};

}

#endif /* ARIADNE5_DipoleSwinger_H */
