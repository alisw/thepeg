// -*- C++ -*-
#ifndef TheP8I_Ropewalk_H
#define TheP8I_Ropewalk_H
//
// This is the declaration of the Ropewalk class.
//

#include "TheP8I/Config/TheP8I.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/ColourSinglet.h"

namespace std {

template <>
struct less<ThePEG::ColourSinglet *> :
    public binary_function<ThePEG::ColourSinglet *,
                           ThePEG::ColourSinglet *, bool> 
{
  /**
   * This is the function called when comparing two pointers to
   * type_info.
   */
  bool operator()(const ThePEG::ColourSinglet * x,
		  const ThePEG::ColourSinglet * y) const {
    if ( x && y ) return x->partons() < y->partons();
    return !y;
  }
};

}

namespace TheP8I {

using namespace ThePEG;

/**
 * Here is the documentation of the Ropewalk class.
 */
class Ropewalk {

public:

  double lambdaBefore;
  vector<double> mspec;
  /**
   * Forward declaration of a Dipole.
   */
  struct Dipole;


  /**
   * Container class storing information of a dipole overlapping with another;
   */
  struct OverlappingDipole {

    /**
     * Standard constructor.
     */
    OverlappingDipole(const Dipole & d,
		      const LorentzRotation R, const Ropewalk * rw);

    /**
     * Check if ovelapping at specific rapidity.
     */
    bool overlap(double y, LorentzPoint b1, Length R0) {
      if ( y < min(ya, yc) || y > max(ya, yc) ) return false;
      LorentzPoint b2 = ba + (bc - ba)*(y - ya)/(yc - ya);
      return (b1 - b2).perp() <= R0;
    }

    /**
     * Pointer to the dipole.
     */
    const Dipole * dipole;

    /**
     * The boosted rapidities.
     */
    double yc, ya;

    /**
     * The propagated and boosted positions.
     */
    LorentzPoint bc, ba;

    /**
     * Relative direction.
     */
    int dir;

    /**
     * Initially estimated fraction of overlap.
     */
    double yfrac;

  };

  /**
   * Container class containing information about overlaps for a dipole.
   * *** ATTENTION *** the meaning of pa and pc is reversed: pa is
   *      always the COLOURED and pc is always the ANTI-COLOURED
   */
  struct Dipole {

    /**
     * Constructor.
     */
    Dipole(tcPPtr p1, tcPPtr p2, ColourSinglet & str)
      : pc(p2), pa(p1), string(&str), n(0), m(0), p(1), q(0), nb(0), broken(false), hadr(false) {}

    /**
     * Return true if this dipole is breakable. Only dipoles between
     * gluons (not belonging to other dipoles which have broken) and
     * which has a mass above 2 GeV can break.
     */
    bool breakable() const {
      return !broken && pc->id() == ParticleID::g && pa->id() == ParticleID::g &&
	!pc->decayed() && !pa->decayed() && s() > 4.0*GeV2;
    }

    /**
     * Calculate the probability that this dipole should break. Taking
     * into account the number of negative steps in the colour space
     * and the fraction of overlapping dipoles which are able to
     * break.
     */
    double breakupProbability() const;

    /**
     * Recalculate overlaps assuming a position a fraction ry from
     * the coloured parton.
     */
    double reinit(double ry, Length R0, Energy m0);

    /**
     * Set and return the effective multiplet.
     */
    void initMultiplet();
    void initNewMultiplet();
    void initWeightMultiplet();

    /**
     * Calculate the multiplicity given pp and qq;
     */
    static double multiplicity(int pp, int qq) {
      return ( pp < 0 || qq < 0 || pp + qq == 0 )? 0.0:
	0.5*(pp + 1)*(qq + 1)*(pp + qq + 2);
    }

    /**
     * Calculate the average kappa-enhancement.
     */
    double kappaEnhancement() const {
      // return 1.0;
      return 0.25*(p + q + 3 - p*q/double(p + q));
    }

    /**
     * Calculate the kappa-enhancement in the first break-up.
     */
    double firstKappa() const {
      return 0.25*(2 + 2*p + q);
    }

    /**
     * Calculate boosted kappa-enhancement.
     */
    double firstKappa(double alpha) const {
      if (alpha > 0)
	return alpha*firstKappa();
      return firstKappa();
    }

    /**
     * Return the squared invariant mass of this dipole.
     */
    Energy2 s() const {
      return (pc->momentum() + pa->momentum()).m2();
    }

    /**
     * Return the effective rapidityspan of this dipole.
     */
    double yspan(Energy m0) const {
      return s() > sqr(m0)? 0.5*log(s()/sqr(m0)): 0.0;
    }

    /**
     * The coloured particle.
     */
    tcPPtr pc;

    /**
     * The anti-coloured particle.
     */
    tcPPtr pa;

    /**
     * The propagated and boosted positions.
     */
    LorentzPoint bc, ba;

    /**
     * The overlapping dipoles with the amount of overlap (negative
     * means anti overlap).
     */
    vector<OverlappingDipole> overlaps;

    /**
     * The string to which this Dipole belongs.
     */
    ColourSinglet * string;

    /**
     * The summed parallel (n) and anti-parallel (m) overlaps from
     * other dipoles.
     */
    int n, m;

    /**
     * The multiplet.
     */
    int p, q;

    /**
     * The number of baryons from junctions
     */
    int nb;

    /**
     * Indicate that this dipole has been broken.
     */
    mutable bool broken;
    mutable bool hadr;
  };


  /**
   * Convenient typedef
   */
  typedef map<ColourSinglet *, vector<Dipole *> > DipoleMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Ropewalk(const vector<ColourSinglet> & singlets, Length width,
	   Energy scale, double jDC, bool throwaway = true, bool rapidityOverlap = true, bool deb = false);

  /**
   * The destructor.
   */
  virtual ~Ropewalk();
  //@}

  /**
   * Return all ColourSinglet objects with their kappa enhancement.
   */
  vector< pair<ColourSinglet,double> > getSinglets(double DeltaY = 0.0) const;

  /**
   * Calculate the average kappa-enhancement for the given colour
   * singlet, optionally only taking into account dipoles the
   * rapidity interval |y| < deltaY. Note that \a cs must point to a
   * object in the stringdipoles map. @return -1 if cs is not found.
   */
  double averageKappaEnhancement(ColourSinglet * cs, double deltaY = 0.0) const {
    DipoleMap::const_iterator it = stringdipoles.find(cs);
    if ( it == stringdipoles.end() ) return -1.0;
    return averageKappaEnhancement(it, deltaY);
  }

  /**
   * Calculate the average kappa-enhancement for the given iterator
   * pointing into the stringdipoles map, optionally only taking
   * into account dipoles the rapidity interval |y| < deltaY. Note
   * that \a cs must point to a object in the stringdipoles
   * map. @return -1 if something went wrong.
   */
  double averageKappaEnhancement(DipoleMap::const_iterator it,
				 double deltaY = 0.0) const {
    return averageKappaEnhancement(it->second, deltaY);
  }

  /**
   * Calculate the average kappa-enhancement the given vector of
   * dipoles, optionally only taking into account dipoles the
   * rapidity interval |y| < deltaY. @return -1 if something went
   * wrong.
   */
  double averageKappaEnhancement(const vector<Dipole *> & dipoles,
				 double deltaY = 0.0) const;

  /**
   * Propagate a parton a finite time inthe lab-system.
   */
  LorentzPoint propagate(const LorentzPoint & b,
			 const LorentzMomentum & p) const;

  /**
  / Get number of juctions
  */
  double getNb();
  
  double getkW();
  /**
  / Get m,n for all strings
  */
  pair<double,double> getMN();
  

  /**
  / Get lambda measure for event
  */
  double lambdaSum();

  /**
   * Calculate the rapidity of a Lorentz vector but limit the
   * transverse mass to be above m0.
   */
  double limitrap(const LorentzMomentum & p) const;

  /**
   * Sort the given vector of dipoles so that the particle \a p is
   * first and the others follows in sequence.  If \a p is a
   * (anti-)triplet, it will be the pc(pa) of the first dipole, and
   * the pc(pa) of the second will be the pa(pc) of the first and so
   * on. If \a p is an octet, the direction will be as if it was a
   * triplet.
   */
  static void sort(vector<Dipole *> & dips, tcPPtr p);

protected:

  /**
   * Extract all dipoles.
   */
  void getDipoles();

  /**
   * Calculate all overlaps and initialize multiplets in all dipoles.
   */
  void calculateOverlaps();

  /**
   * Break dipoles (and strings)
   */
  void doBreakups();

  /**
   * Return the current step.
   */
  Step & step() const;

  /**
   * Make a copy of a ColourSinglet making sure all partons are final.
   */
  static ColourSinglet * cloneToFinal(const ColourSinglet & cs);

public:

  DipoleMap::iterator begin() {
    return stringdipoles.begin();
  }

  DipoleMap::iterator end() {
    return stringdipoles.end();
  }

  /**
   * Exception class.
   */
  struct ColourException: public Exception {};

private:

  /**
   * The assumed radius of a string.
   */
  Length R0;

  /**
   * The assumed mass for calculating rapidities
   */
  Energy m0;

  double junctionDiquarkProb;

  /**
   * The vector of all Dipoles
   */
  vector<Dipole> dipoles;

  /**
   * All Dipoles in all string
   */
  DipoleMap stringdipoles;


  /**
   * Do rapidity overlap=
   */
  bool rapidityOverlap;
  /**
   * Debug flag.
   */
  bool debug;

  public:

  mutable double strl0;
  mutable double strl;
  mutable double avh;
  mutable double avp;
  mutable double avq;
  

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Ropewalk & operator=(const Ropewalk &);

  double nb;
};

}

#endif /* TheP8I_Ropewalk_H */
