// Header file for the Rndm, Hist, Vec4 and Particle classes,
// and for some convenient global functions.

#ifndef Basics_H
#define Basics_H

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include "Pythia7/Config/Pythia7.h"

namespace Pythia7 {

/**
 * The underlying classes in the Pythia7 parton shower handler do not
 * directly use the ThePEG framework. All these classes are included
 * in the Shower namespace.
 */
namespace Shower {

using namespace std;

//**************************************************************************

// Some simple global functions that are very convenient.

// Minimum and maximum for long integers - should have been part of compiler?
// inline long min(long i, long j) {return (i < j) ? i : j;}
// inline long max(long i, long j) {return (i > j) ? i : j;}

// Ensure that roundoff does not give square root of negative number.

/**
 * Return the square root. If the argument is negative, zero will be
 * returned.
 */
inline double sqrtpos(double x) {if (x < 0.) return 0.; else return sqrt(x);}

/**
 * Return the square root. If the argument \a x is negative the
 * negative of the square root of the absolute, \f$-\sqrt{-x}\f$, is
 * returned.
 */
inline double sqrtsgn(double x) {if (x < 0.) return -sqrt(-x); 
  else return sqrt(x);}

/**
 * The Kallens lambda function. \f$\lambda_K(a,b,c)=(a-b-c)^2-4bc\f$.
 */
inline double lambdaK(double a, double b, double c) {
  return max(0., pow(a - b - c, 2) - 4. * b * c);
}

/**
 * Square root of the Kallens lambda
 * function. \f$\sqrt{\lambda_K(a,b,c)}=\sqrt{(a-b-c)^2-4bc}\f$.
 */
inline double lambdaKRoot(double a, double b, double c) {
  return sqrt(lambdaK(a, b, c));
}

/**
 * Magnitude of three-momentum. Given the mass squared, \a a, of a
 * decaying particle return the magnitude of the three-momentum of the
 * two daughters with squared masses \a b and \a c in the rest system
 * of the decay.
 */
inline double pAbsLambdaK(double a, double b, double c) {
  return 0.5*sqrt(lambdaK(a, b, c)/a);
}

/**
 * The random number generator class used internally in the Pythia7
 * Shower routines. This should be replaced with the ThePEG
 * RandomGenerator. Uses only static methods.
 */
class Rndm {

public:

  /** Default constructor. */
  Rndm() {;} 
  /**
   * Standard constructor giving a seed. Note that although this looks
   * like a new Rndm object is constructed with a separate random
   * number sequence, the underlying static generator is
   * reinitialized.
   */
  Rndm(long inseed) {init(inseed);}

  /**
   * Initialize the random generator optionally giving a seed.
   */
  static void init(long = 0);

  /**
   * Return a number uniformly distributed between 0 and 1.
   */
  static double flat() ;

  /**
   * Return a random number distributed according to \f$\exp(x)\f$.
   */
  static double exp() { return  -log(flat()) ;} 

  /**
   * Return a random number distributed according to \f$x\exp(x)\f$.
   */
  static double xexp() { return  -log(flat() * flat()) ;} 

  /**
   * Return a random number distributed according to
   * \f$\exp(-x^2/2)\f$.
   */
  static double gauss() ;

  /**
   * Return an index of the given vector according to the relative
   * (positive) probabilities therein.
   */
  static long pick(const vector<double>&) ; 

private:

  /**
   * Has the generator been initialized?
   */
  static bool initrndm;

  /**
   * Has a gaussian number been
   * saved?
   */
  static bool savgauss;

  /**
   * Internal index.
   */
  static long i97;

  /**
   * Internal index.
   */
  static long j97;

  /**
   * Internal variable.
   */
  static double u[97];

  /**
   * Internal variable.
   */
  static double c;

  /**
   * Internal variable.
   */
  static double cd;

  /**
   * Internal variable.
   */
  static double cm;

  /**
   * Internal variable.
   */
  static double save;
};

//**************************************************************************

/**
 * Simple histogram class to be used internally in the Pythia7 Shower
 * classes for testing purposes.
 */
class Hist{

public:

  /** The maximum number of lines a histogram can use at output. */
  static long NLINES;

  /** Tolerance in deviation of xmin and xmax between two histograms. */
  static double TOLERANCE;

  /** Small number to avoid division by zero. */
  static double TINY;

  /**
   * Default constructor.
   */
  Hist() {}

  /**
   * Standard constructor.
   * @param titlein the title of the histogram.
   * @param nbinin the number of bins.
   * @param xminin the lower edge of the histogram.
   * @param xmaxin the upper edge of the histogram
   */
  Hist(string titlein, long nbinin = 100, double xminin = 0., 
    double xmaxin = 1.) {
    book(titlein, nbinin, xminin, xmaxin);
  }

  /**
   * Copy constructor.
   */ 
  Hist(const Hist& h) {
    title = h.title; nbin = h.nbin; xmin = h.xmin; xmax = h.xmax; 
    dx = h.dx; nfill = h.nfill; under = h.under; inside = h.inside; 
    over = h.over; res = h.res;
  }
  
  /**
   * Copy constructor giving a new title.
   */    
  Hist(string titlein, const Hist& h) {
    title = titlein; nbin = h.nbin; xmin = h.xmin; xmax = h.xmax; 
    dx = h.dx; nfill = h.nfill; under = h.under; inside = h.inside; 
    over = h.over;res = h.res; 
  }

  /**
   * Assignment operator.
   */
  Hist& operator=(const Hist& h) {
    if(this != &h) {
      nbin = h.nbin; xmin = h.xmin; xmax = h.xmax; 
      dx = h.dx; nfill = h.nfill; under = h.under; inside = h.inside; 
      over = h.over; res = h.res;
      return *this;
    }
  }    

  /**
   * (Re-) Book this histogram. (Same parameters as in the Standard
   * constructor)
   */
  void book(string = "  ", long = 100, double = 0., double = 1.);

  /**
   * Set a new name.
   */
  void name(string = "  ");  

  /**
   * Reset this histogram.
   */
  void null() ; 

  /**
   * Fill histogram.
   */
  void fill(double, double = 1.);

  /**
   * Print histogram to the given stream. Print a table of x and y
   * values suitable to be used by programs such as gnuplot.
   */
  void table(ostream& = cout) const ;

  /**
   * Check if the given histogram is compatible with this.
   */
  bool sameSize(const Hist&) const ;

  /** Add histogram. */
  Hist& operator+=(const Hist&) ; 
  /** Subtract histogram. */
  Hist& operator-=(const Hist&) ;
  /** Multiply by histogram (bin-by-bin). */
  Hist& operator*=(const Hist&) ; 
  /** Divide by histogram (bin-by-bin). */
  Hist& operator/=(const Hist&) ;
  /** Add number to each bin. */
  Hist& operator+=(double) ; 
  /** Subtract number from each bin. */
  Hist& operator-=(double) ; 
  /** Multiply by number in each bin. */
  Hist& operator*=(double) ; 
  /** Divide by number in each bin. */
  Hist& operator/=(double) ; 

  /** Add number and histogram. */
  friend Hist operator+(double, const Hist&);
  /** Add number and histogram. */
  friend Hist operator+(const Hist&, double);
  /** Add two histograms. */
  friend Hist operator+(const Hist&, const Hist&);
  /** Subtraction between number and histogram. */
  friend Hist operator-(double, const Hist&);
  /** Subtraction between number and histogram. */
  friend Hist operator-(const Hist&, double);
  /** Subtraction between two histograms. */
  friend Hist operator-(const Hist&, const Hist&);
  /** Multiply number with histogram (bin-by-bin). */
  friend Hist operator*(double, const Hist&);
  /** Multiply number with histogram (bin-by-bin). */
  friend Hist operator*(const Hist&, double);
  /** Multiply two histograms (bin-by-bin). */
  friend Hist operator*(const Hist&, const Hist&);
  /** Divide number with histogram (bin-by-bin). */
  friend Hist operator/(double, const Hist&);
  /** Divide histogram with number (bin-by-bin). */
  friend Hist operator/(const Hist&, double);
  /** Divide two histograms (bin-by-bin). */
  friend Hist operator/(const Hist&, const Hist&);

  /** Standard output operator. */
  friend ostream& operator<<(ostream&, const Hist&) ;

  /** Non-member function of table(ostream&). */
  friend void table(const vector<Hist>&, ostream& = cout) ;

private:

  /** The title. */
  string title;
  /** The number of bins. */
  long nbin;
  /** The number of times fill() has been called. */
  long nfill;
  /** the lower cut. */
  double xmin; 
  /** the upper cut. */
  double xmax; 
  /** the bin width. */
  double dx; 
  /** the underflow bin. */
  double under; 
  /** the sum of all bins. */
  double inside; 
  /** the overflow bin. */
  double over; 
  /** The bins. */
  vector<double> res;

};

//**************************************************************************

// Forward reference to RotBstMatrix class.
class RotBstMatrix;

//**************************************************************************


/**
 * Vec4 class.
 * This class implements four-vectors, in energy-momentum space.
 * (But could equally well be used to hold space-time four-vectors.)
 * Used by the internal Pythia7 Shower classes.
 */
class Vec4 {

public:
  /**
   * Constructors.
   */
  Vec4(double xin = 0., double yin = 0., double zin = 0., double tin = 0.)
    : xx(xin), yy(yin), zz(zin), tt(tin) {}
  /** NOT DOCUMENTED */
  Vec4(const Vec4& v)
    : xx(v.xx), yy(v.yy), zz(v.zz), tt(v.tt) {}
  /** NOT DOCUMENTED */
  Vec4& operator=(const Vec4& v) {
    if (this != &v) {xx = v.xx; yy = v.yy; zz = v.zz; tt = v.tt; }
    return *this; }

      
  /**
   * Member functions for input.
   */
  void p(double xin, double yin, double zin, double tin) 
    {xx = xin; yy = yin; zz = zin; tt = tin;}
  /** NOT DOCUMENTED */
  void px(double xin) {xx = xin;}
  /** NOT DOCUMENTED */
  void py(double yin) {yy = yin;}
  /** NOT DOCUMENTED */
  void pz(double zin) {zz = zin;}
  /** NOT DOCUMENTED */
  void e(double tin) {tt = tin;}


  /**
   * Member functions for output.
   */
  double px() const {return xx;}
  /** NOT DOCUMENTED */
  double py() const {return yy;}
  /** NOT DOCUMENTED */
  double pz() const {return zz;}
  /** NOT DOCUMENTED */
  double e() const {return tt;}
  /** NOT DOCUMENTED */
  double m2calc() const {return tt*tt - xx*xx - yy*yy - zz*zz;}
  /** NOT DOCUMENTED */
  double mcalc() const {return sqrtsgn(tt*tt - xx*xx - yy*yy - zz*zz);}
  /** NOT DOCUMENTED */
  double pT() const {return sqrt(xx*xx + yy*yy);}
  /** NOT DOCUMENTED */
  double pT2() const {return xx*xx + yy*yy;}
  /** NOT DOCUMENTED */
  double pAbs() const {return sqrt(xx*xx + yy*yy + zz*zz);}
  /** NOT DOCUMENTED */
  double p2() const {return xx*xx + yy*yy + zz*zz;}
  /** NOT DOCUMENTED */
  double theta() const {return atan2(sqrt(xx*xx + yy*yy), zz);}
  /** NOT DOCUMENTED */
  double phi() const {return atan2(yy,xx);}


  /**
   * Member functions that perform operations.
   */
  void rescalep(double fac) {xx = fac * xx; yy = fac * yy; zz = fac * zz;}
  /** NOT DOCUMENTED */
  void flip(){xx = -xx; yy = -yy; zz = -zz;}
  /** NOT DOCUMENTED */
  void rot(double, double); 
  /** NOT DOCUMENTED */
  void rotaxis(double, double, double, double); 
  /** NOT DOCUMENTED */
  void rotaxis(double, const Vec4&);
  /**
   * Member functions that perform operations.
   */
  void bst(double, double, double); 
  /** NOT DOCUMENTED */
  void bst(const Vec4&); 
  /** NOT DOCUMENTED */
  void rotbst(const RotBstMatrix&); 


  /**
   * Operator overloading with member functions
   */
  Vec4& operator+=(const Vec4& v) {xx += v.xx; yy += v.yy; zz += v.zz; 
    tt += v.tt; return *this;}
  /** NOT DOCUMENTED */
  Vec4& operator-=(const Vec4& v) {xx -= v.xx; yy -= v.yy; zz -= v.zz; 
    tt -= v.tt; return *this;}
  /** NOT DOCUMENTED */
  Vec4& operator*=(double f) {xx *= f; yy *= f; zz *= f; 
    tt *= f; return *this;}
  /** NOT DOCUMENTED */
  Vec4& operator/=(double f) {xx /= f; yy /= f; zz /= f; 
    tt /= f; return *this;}


  /**
   * Operator overloading with friends
   */
  friend Vec4 operator+(const Vec4&, const Vec4&);
  /**
   * Operator overloading with friends
   */
  friend Vec4 operator-(const Vec4&, const Vec4&);
  /**
   * Operator overloading with friends
   */
  friend Vec4 operator*(double, const Vec4&);
  /**
   * Operator overloading with friends
   */
  friend Vec4 operator*(const Vec4&, double);
  /**
   * Operator overloading with friends
   */
  friend Vec4 operator/(const Vec4&, double);
  /**
   * Operator overloading with friends
   */
  friend double operator*(const Vec4&, const Vec4&);

  /**
   * Scalar product of 3-vector parts.
   */
  friend double dot3(const Vec4&, const Vec4&);

  /**
   * Cross product of 3-vector parts.
   */
  friend Vec4 cross3(const Vec4&, const Vec4&);

  /**
   * Cosine of the polar angle between \a v1 and \a v2.
   */
  friend double costheta(const Vec4 & v1, const Vec4 & v2);

  /**
   * The polar angle between \a v1 and \a v2.
   */
  friend double theta(const Vec4& v1, const Vec4& v2) {
    return acos(costheta(v1, v2)); } 

  /**
   * Cosine of the azimuthal angle between \a v1 and \a v2 around the
   * \a n axis.
   */
  friend double cosphi(const Vec4 & v1, const Vec4 & v2, const Vec4 & n);

  /**
   * The azimuthal angle between \a v1 and \a v2 around the \a n axis.
   */
  friend double phi(const Vec4 & v1, const Vec4 & v2, const Vec4 & n) {
    return acos(cosphi(v1, v2, n)); } 


  /**
   * Print a four-vector
   */
  friend ostream& operator<<(ostream&, const Vec4&) ;

private:

  /** The x component. */
  double xx;

  /** The y component. */
  double yy;

  /** The z component. */
  double zz;

  /** The time component. */
  double tt;

};

// Implementation of operator overloading with friends.

/** NOT DOCUMENTED */
inline Vec4 operator+(const Vec4& v1, const Vec4& v2) 
  {Vec4 v = v1 ; return v += v2;}

/** NOT DOCUMENTED */
inline Vec4 operator-(const Vec4& v1, const Vec4& v2) 
  {Vec4 v = v1 ; return v -= v2;}

/** NOT DOCUMENTED */
inline Vec4 operator*(double f, const Vec4& v1) 
  {Vec4 v = v1; return v *= f;}

/** NOT DOCUMENTED */
inline Vec4 operator*(const Vec4& v1, double f) 
  {Vec4 v = v1; return v *= f;}

/** NOT DOCUMENTED */
inline Vec4 operator/(const Vec4& v1, double f) 
  {Vec4 v = v1; return v /= f;}

/** NOT DOCUMENTED */
inline double operator*(const Vec4& v1, const Vec4& v2)
  {return v1.tt*v2.tt - v1.xx*v2.xx - v1.yy*v2.yy - v1.zz*v2.zz;}  

//**************************************************************************


/**
 * RotBstMatrix class.
 * This class implements 4 * 4 matrices that encode an arbitrary combination
 * of rotations and boosts, that can be applied to Vec4 four-vectors.
 * Used by the internal Pythia7 Shower classes.
 */
class RotBstMatrix {

public:
  /**
   * Constructors.
   */
  RotBstMatrix() {for (long i = 0; i < 4; ++i) { for (long j = 0; j < 4; ++j) {
    M[i][j] = (i==j) ? 1. : 0.; } } } 
  /** NOT DOCUMENTED */
  RotBstMatrix(const RotBstMatrix& Min) {
    for (long i = 0; i < 4; ++i) { for (long j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } } }
  /** NOT DOCUMENTED */
  RotBstMatrix& operator=(const RotBstMatrix& Min) {if (this != &Min) {
    for (long i = 0; i < 4; ++i) { for (long j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } }} return *this; }


  /**
   * Member functions.
   */
  void rot(double = 0., double = 0.);
  /**
   * Member functions.
   */
  void rot(const Vec4& p);
  /**
   * Member functions.
   */
  void bst(double = 0., double = 0., double = 0.);
  /**
   * Member functions.
   */
  void bst(const Vec4&);
  /**
   * Member functions.
   */
  void bstback(const Vec4&);
  /**
   * Member functions.
   */
  void bst(const Vec4&, const Vec4&);
  /**
   * Member functions.
   */
  void rotbst(const RotBstMatrix&);
  /**
   * Member functions.
   */
  void invert();
  /**
   * Member functions.
   */
  void toCMframe(const Vec4&, const Vec4&);
  /**
   * Member functions.
   */
  void fromCMframe(const Vec4&, const Vec4&);
  /**
   * Member functions.
   */
  void reset();

  /**
   * Print a transformation matrix.
   */
  friend ostream& operator<<(ostream&, const RotBstMatrix&) ;

  /**
   * Private members to be accessible from Vec4. 
   */
  friend class Vec4;

private:
  /** NOT DOCUMENTED */
  double M[4][4];

};

//**************************************************************************


/**
 * Particle class.
 * This class holds info on a particle in general.
 * Used by the internal Pythia7 Shower classes.
 */
class Particle {

public:
  /**
   * Constructors.
   */
  Particle() {idp = 0; statusp = 0; mother1p = -1; mother2p = -1; prevp = -1; 
    colp = 0; anticolp = 0; pp = Vec4(0.,0.,0.,0.); mp = 0.; scalep = 0.;}
  /** NOT DOCUMENTED */
  Particle(long idin, long statusin = 0, long mother1in = -1, 
    long mother2in = -1, long colin = 0, long anticolin = 0, 
    Vec4 pin = Vec4(0.,0.,0.,0.), double min = 0., double scalein = 0.) { 
    idp = idin; statusp = statusin; mother1p = mother1in; 
    mother2p = mother2in; prevp = -1; colp = colin; anticolp = anticolin; 
    pp = pin; mp = min; scalep = scalein;}  
  /** NOT DOCUMENTED */
  Particle(long idin, long statusin = 0, long mother1in = -1, 
    long mother2in = -1, long colin = 0, long anticolin = 0, 
    double pxin = 0., double pyin = 0., double pzin = 0.,
    double ein = 0., double min = 0., double scalein = 0.) { 
    idp = idin; statusp = statusin; mother1p = mother1in; 
    mother2p = mother2in; prevp = -1; colp = colin; anticolp = anticolin; 
    pp = Vec4(pxin, pyin, pzin, ein); mp = min; scalep = scalein;}  
  /** NOT DOCUMENTED */
  Particle(const Particle& pt) {
    idp = pt.idp; statusp = pt.statusp; mother1p = pt.mother1p; 
    mother2p = pt.mother2p; prevp = pt.prevp;
  /**
   * Constructors.
   */
    colp = pt.colp; anticolp = pt.anticolp; 
    pp = pt.pp; mp = pt.mp; scalep = pt.scalep;} 
  /** NOT DOCUMENTED */
  Particle& operator=(const Particle& pt) {if (this != &pt) {
    idp = pt.idp; statusp = pt.statusp; mother1p = pt.mother1p; 
    mother2p = pt.mother2p; prevp = pt.prevp;
    colp = pt.colp; anticolp = pt.anticolp; 
    pp = pt.pp; mp = pt.mp; scalep = pt.scalep;} return *this; } 

      
  /**
   * Member functions for input.
   */
  void id(long idin) {idp = idin;}
  /** NOT DOCUMENTED */
  void status(long statusin) {statusp = statusin;}
  /** NOT DOCUMENTED */
  void addstatus(long change) {statusp += change;}
  /** NOT DOCUMENTED */
  void mother1(long mother1in) {mother1p = mother1in;}
  /** NOT DOCUMENTED */
  void mother2(long mother2in) {mother2p = mother2in;}
  /** NOT DOCUMENTED */
  void mothers(long mother1in = -1, long mother2in = -1) 
    {mother1p = mother1in; mother2p = mother2in;}
  /** NOT DOCUMENTED */
  void prev(long previn) {prevp = previn;}
  /** NOT DOCUMENTED */
  void col(long colin) {colp = colin;}
  /** NOT DOCUMENTED */
  void anticol(long anticolin) {anticolp = anticolin;}
  /** NOT DOCUMENTED */
  void cols(long colin = 0,long anticolin = 0) 
    {colp = colin; anticolp = anticolin;}  
  /** NOT DOCUMENTED */
  void p(Vec4 pin) {pp = pin;}
  /** NOT DOCUMENTED */
  void p(double pxin, double pyin, double pzin, double ein) 
    {pp.p(pxin, pyin, pzin, ein);}
  /** NOT DOCUMENTED */
  void px(double pxin) {pp.px(pxin);}
  /** NOT DOCUMENTED */
  void py(double pyin) {pp.py(pyin);}
  /** NOT DOCUMENTED */
  void pz(double pzin) {pp.pz(pzin);}
  /** NOT DOCUMENTED */
  void e(double ein) {pp.e(ein);}
  /** NOT DOCUMENTED */
  void m(double min) {mp = min;}
  /** NOT DOCUMENTED */
  void scale(double scalein) {scalep = scalein;}


  /**
   * Member functions for output.
   */
  long id() const {return idp;}
  /** NOT DOCUMENTED */
  long status() const {return statusp;}
  /** NOT DOCUMENTED */
  long mother1() const {return mother1p;}
  /** NOT DOCUMENTED */
  long mother2() const {return mother2p;}
  /** NOT DOCUMENTED */
  long prev() const {return prevp;}
  /** NOT DOCUMENTED */
  long col() const {return colp;}
  /** NOT DOCUMENTED */
  long anticol() const {return anticolp;}
  /** NOT DOCUMENTED */
  Vec4 p() const {return pp;}
  /** NOT DOCUMENTED */
  double px() const {return pp.px();}
  /** NOT DOCUMENTED */
  double py() const {return pp.py();}
  /** NOT DOCUMENTED */
  double pz() const {return pp.pz();}
  /** NOT DOCUMENTED */
  double e() const {return pp.e();}
  /** NOT DOCUMENTED */
  double m() const {return mp;}
  /** NOT DOCUMENTED */
  double scale() const {return scalep;}
  /** NOT DOCUMENTED */
  double m2() const {return mp*mp;}
  /** NOT DOCUMENTED */
  double mcalc() const {return pp.mcalc();}
  /** NOT DOCUMENTED */
  double m2calc() const {return pp.m2calc();}
  /** NOT DOCUMENTED */
  double pT() const {return pp.pT();}
  /** NOT DOCUMENTED */
  double pT2() const {return pp.pT2();}
  /** NOT DOCUMENTED */
  double mT() const {return sqrt(mp*mp + pp.pT2());}
  /** NOT DOCUMENTED */
  double mT2() const {return mp*mp + pp.pT2();}
  /** NOT DOCUMENTED */
  double pAbs() const {return pp.pAbs();}
  /** NOT DOCUMENTED */
  double p2() const {return pp.p2();}
  /** NOT DOCUMENTED */
  double theta() const {return pp.theta();}
  /** NOT DOCUMENTED */
  double phi() const {return pp.phi();}


  /**
   * Member functions that perform operations.
   */
  void rescalep(double fac) {pp.rescalep(fac);}
  /** NOT DOCUMENTED */
  void rot(double theta, double phi) {pp.rot(theta, phi);} 
  /** NOT DOCUMENTED */
  void bst(double betaX, double betaY, double betaZ) 
    {pp.bst(betaX, betaY, betaZ);}
  /** NOT DOCUMENTED */
  void bst(const Vec4& vec) {pp.bst(vec);}
  /** NOT DOCUMENTED */
  void rotbst(const RotBstMatrix& M) {pp.rotbst(M);} 


  /**
   * Print a particle
   */
  friend ostream& operator<<(ostream&, const Particle&) ;

private:
  /** NOT DOCUMENTED */
  long idp;
  /** NOT DOCUMENTED */
  long statusp;
  /** NOT DOCUMENTED */
  long mother1p;
  /** NOT DOCUMENTED */
  long mother2p;
  /** NOT DOCUMENTED */
  long prevp;
  /** NOT DOCUMENTED */
  long colp;
  /** NOT DOCUMENTED */
  long anticolp;
  /** NOT DOCUMENTED */
  Vec4 pp;
  /** NOT DOCUMENTED */
  double mp;
  /** NOT DOCUMENTED */
  double scalep;
};

//**************************************************************************

/** Function to give particle masses. */
double Mass(long);

/** Function to give spin of particle. */
long iSpin(long);

/** Function to give charge of particle. */
long iCharge(long);

/** Function to give colour of particle. */
long iColour(long);

//**************************************************************************
}
}

#endif // Basics_H
