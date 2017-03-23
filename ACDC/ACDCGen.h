// -*- C++ -*-
//
// ACDCGen.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ACDCGen_H
#define ACDCGen_H

#include <algorithm>
#include "ACDCGenConfig.h"
#include "ACDCTraits.h"
#include "ACDCGenCell.h"
#include "ThePEG/Utilities/Exception.h"

namespace ACDCGenerator {

/**
 * ACDCGen is a general class for sampling multi-dimensional
 * functions. ACDCGen can sample several functions simultaneously,
 * selecting different functions according to the relative
 * probabilities determined by their total integrals. The functions
 * are sampled on a unit hypercube. Function object of any class can
 * be used as long as the ACDCFncTraits class is specialized
 * correctly. ACDCFncTraits can also be used to rescale values in the
 * unit hypercube to any desired range. ACDCGen needs a random number
 * generator. Again, random number generators of any class can be used
 * as long as the ACDCRandomTraits class is specialized correctly.
 *
 * To give an unweighted samlpe ACDCGen uses a compensating
 * algorithm. Before the production sampling begins, the functions are
 * sampled randomly in the hypercube a user-defined number of times to
 * find an approximate maxumum value. The hypercube is then divided
 * into cells each of which have an approximate maximum value of the
 * function, to enable efficient sampling. The maxima are only
 * approximate though and, if a function value is found above the
 * maximum in a cell the ACDCGen will go into a compensating mode. The
 * cell is then first subdivided further and in the following this
 * cell will be over-sampled to compensate for that fact that it was
 * under-sampled before. In this way the probability of obtaining a
 * biased sample is reduced. Also rather functions with large peaks
 * are then sampled rather efficiently. Functions with narrow peaks
 * should, however, be avoided since there is no guarantee that the
 * peack is actually hit.
 */
template <typename Rnd, typename FncPtr>
class ACDCGen {

public:

  /** Template argument typedef. */
  typedef Rnd RndType;
  /** Template argument typedef. */
  typedef ACDCRandomTraits<RndType> RndTraits;
  /** Template argument typedef. */
  typedef FncPtr FncPtrType;
  /** A vector of cells. */
  typedef vector<ACDCGenCell*> CellVector;
  /** A vector of function objects. */
  typedef vector<FncPtrType> FncVector;
  /** A vector of integers. */
  typedef vector<DimType> DimVector;
  /** The size type of the vectors used. */
  typedef DimVector::size_type size_type;
  /** Template argument typedef. */
  typedef ACDCFncTraits<FncPtrType> FncTraits;

public:

  /**
   * Standard constructor requiring a random generator object to be
   * used.
   */
  inline ACDCGen(Rnd * r);

  /**
   * Default Constructor.
   */
  inline ACDCGen();

  /**
   * Destructor.
   */
  inline ~ACDCGen();

  /**
   * Add a function of a given dimension, \a dim, according to which
   * points will be generated. Note that each function, \a f, added
   * like this will have its own tree of cells. The \a maxrat argument
   * determines the lowest ratio of values allowed between the cell
   * with lowest and highest value. If negative it is given by 1/nTry().
   */
  inline bool addFunction(DimType dim, FncPtrType f, double maxrat = -1.0);

  /**
   * Remove all added functions and reset the generator;
   */
  inline void clear();

public:

  /**
   * Generate a point, choosing between the different functions
   * specified. The chosen function is returned, while the generated
   * point is obtained by the function lastPoint().
   */
  inline FncPtrType generate();

  /**
   * Reject the last generated point. Only used in the evaluation of
   * the total integral.
   */
  inline void reject();

  /**
   * Return the last generated point.
   * @return a vector of doubles, each in the interval ]0,1[.
   */
  inline const DVector & lastPoint() const;

  /**
   * Return the value of the last chosen function in the last point.
   */
  inline double lastF() const;

  /**
   * Return the function chosen for the last generated point.
   */
  inline FncPtrType lastFunction() const;

  /**
   * return the index of the function chosen for the last generated
   * point.
   */
  inline size_type last() const;

public:

  /** @name Functions influencing the efficiency of the generation. */
  //@{
  /**
   * Set the minimum cell size considered for this generation. The
   * default is the machine limit for double precision times a
   * hundred.
   */
  inline void eps(double newEps);

  /**
   * Set the safety margin used to multiply the highest found function
   * value in a cell when setting its overestimated value. (Default is
   * 1.1.)
   */
  inline void margin(double newMargin);

  /**
   * Set the number of points (with non-zero function value) used to
   * initialize the tree of cells to use in the generation for each
   * function.
   */
  inline void nTry(size_type newNTry);

  /**
   * Set the maximum number of attempts to generate a phase space
   * point, or to find non-zero points in the initialization.
   */
  inline void maxTry(long);
  //@}

public:

  /** @name Information about the current generation. */
  //@{
  /**
   * Return the current Monte Carlo estimate of the integral of the
   * specified function (or all functions if NULL) over the unit volume. 
   */
  inline double integral(FncPtrType f = FncPtrType()) const;

  /**
   * Return the error on the current Monte Carlo estimate of the
   * integral of the specified function (or all functions if NULL)
   * over the unit volume.
   */
  inline double integralErr(FncPtrType f = FncPtrType()) const;

  /**
   * The number of accepted points so far.
   */
  inline long n() const;

  /**
   * The number of calls to generate() so far. Note that the number of
   * calls to the specified functions may be larger. It is up to the
   * user to keep track of those.
   */
  inline long N() const;

  /**
   * The ratio of the number of accepted and number of tried points
   * n()/N();
   */
  inline double efficiency() const;

  /**
   * Return the number of active cells created so far.
   */
  inline int nBins() const;

  /**
   * Return the maximum depth of any tree of cells used.
   */
  inline int depth() const;

  /**
   * Return the current overestimation of the full integral of all
   * specified functions over the unit volume.
   */
  inline double maxInt() const;
  //@}

  /** @name Access to member variables. */
  //@{
  /**
   * The minimum cell size considered for this generation.
   */
  inline double eps() const;

  /**
   * The safety margin used to multiply the highest found function
   * value in a cell when setting its overestimated value.
   */
  inline double margin() const;

  /**
   * The number of points used to initialize the tree of cells to use
   * in the generation.
   */
  inline size_type nTry() const;

  /**
   * The maximum number of attempts to generate a phase space point,
   * or to find non-zero points in the initialization.
   */
  inline long maxTry() const;

  /**
   * Returns true if generating random numbers are so cheap that a new
   * one can be thrown everytime a sub-cell is chosen. Otherwise
   * random numbers used for this will be reused.
   */
  inline bool cheapRandom() const;

  /**
   * The number of functions used.
   */
  inline size_type size() const;

  /**
   * Returns true if the generator is currently in a state of
   * compensating an erroneous overestimation of one of the specified
   * functions. If so, the integral and the generated points are not
   * statistically correct.
   */
  inline bool compensating();

  /**
   * Return an estimate of how many points need to be sampled before
   * the generator finishes compensating.
   */
  inline long compleft() const;

  /**
   * Return a vector with information about all cells.
   */
  vector<ACDCGenCellInfo> extractCellInfo() const;
  //@}

public:

  /** @name Functions related to the random number generator. */
  //@{
  /**
   * Set to true if generating random numbers are so cheap that a new
   * one can be thrown everytime a sub-cell is chosen. Otherwise
   * random numbers used for this will be reused.
   */
  inline void cheapRandom(bool b);

  /**
   * Set a new random number generator.
   */
  inline void setRnd(Rnd * r);

  /**
   * Double precision number in the interval ]0,1[.
   */
  inline double rnd() const;

  /**
   * Double precision number in the interval ]lo,up[.
   */
  inline double rnd(double lo, double up) const;

  /**
   * Fill the r vector with doubles r[i] in the interval ]lo[i],up[i][.
   */
  inline void rnd(const DVector & lo, const DVector & up, DVector & r)const;

  /**
   * Fill the D first elements in the r vector with doubles in the
   * interval ]0,1[.
   */
  inline void rnd(DimType D, DVector & r) const;

  /**
   * Integer in the interval [0,x[
   */
  inline long rndInt(long x) const;
  //@}

public:

  /**
   * This function is to be used in ThePEG for output to
   * a persistent stream and will not work properly for normal
   * ostreams.
   */
  template <typename POStream>
  void output(POStream &) const;

  /**
   * This function is to be used in ThePEG for input from a persistent
   * stream and will not work properly for normal istreams.
   */
  template <typename PIStream>
  void input(PIStream &);

private:

  /**
   * Calculate the overestimated integral for all functions.
   */
  inline double doMaxInt();

  /**
   * Return the vector of functions.
   */
  inline const FncVector & functions() const;

  /**
   * Return the i'th function.
   */
  inline FncPtrType function(size_type i) const;

  /**
   * Return a vector with the dimensions of all functions.
   */
  inline const DimVector & dimensions() const;

  /**
   * Return the dimension of the i'th function.
   */
  inline DimType dimension(size_type i) const;

  /**
   * Return the dimension of the function chosen for the last
   * generated point.
   */
  inline DimType lastDimension() const;

  /**
   * Return the roots of all cell trees.
   */
  inline const CellVector & cells() const;

  /**
   * Return the root cell for the i'th function.
   */
  inline ACDCGenCell * cell(size_type i) const;

  /**
   * Return the root cell for the function chosen for the last
   * generated point.
   */
  inline ACDCGenCell * lastPrimary() const;

  /**
   * Return a vector with the incremental sum of overestimated
   * integrals for each function.
   */
  inline const DVector & sumMaxInts() const;

  /**
   * Return the cell chosen for the last generated point.
   */
  inline ACDCGenCell * lastCell() const;


  /**
   * Choose a function according to its overestimated integral and
   * choose a cell to generate a point in.
   */
  inline void chooseCell(DVector & lo, DVector & up);

  /**
   * Start the compensation procedure for the last chosen cell when a
   * function velue has been found which exceeds the previous
   * overestimation.
   */
  inline void  compensate(const DVector & lo, const DVector & up);

private:

  /**
   * The random number generator to be used for this Generator.
   */
  RndType * theRnd;
  
  /**
   * The number of accepted points (weight > 0) so far.
   */
  long theNAcc;

  /**
   * The number of attempted points so far.
   */
  long theN;

  /**
   * The number of attempts per function so far.
   */
  vector<long> theNI;

  /**
   * The summed weights per function so far.
   */
  DVector theSumW;

  /**
   * The summed squared weights per function so far.
   */
  DVector theSumW2;

  /**
   * The smallest possible division allowed.
   */
  double theEps;

  /**
   * The factor controlling the loss of efficiency when compensating.
   */
  double theMargin;

  /**
   * The number of points to use to find initial average.
   */
  size_type theNTry;

  /**
   * The maximum number of attempts to generate a phase space point,
   * or to find non-zero points in the initialization.
   */
  long theMaxTry;

  /**
   * True if generating random numbers are so cheap that a new one can
   * be thrown everytime a sub-cell is chosen. Otherwise random
   * numbers used for this will be reused.
   */
  bool useCheapRandom;

  /**
   * A vector of functions.
   */
  FncVector theFunctions;

  /**
   * The dimensions of the functions in theFunctions.
   */
  DimVector theDimensions;

  /**
   * The root of the cell tree for the functions in theFunctions.
   */
  CellVector thePrimaryCells;

  /**
   * The accumulated sum of overestimated integrals of the functions
   * in theFunctions.
   */
  DVector theSumMaxInts;

  /**
   * The last index chosen
   */
  size_type theLast;

  /**
   * The last cell chosen.
   */
  ACDCGenCell * theLastCell;

  /**
   * The last point generated.
   */
  DVector theLastPoint;

  /**
   * The function value of the last point.
   */
  double theLastF;

  /**
   * A helper struct representing a level of compensation.
   */
  struct Level {

    /**
     * The number of attempts after which this level disapprears.
     */
    long lastN;

    /**
     * The previous max value in the Cell to compensate.
     */
    double g;

    /**
     * The cell which is being compensated.
     */
    ACDCGenCell * cell;

    /**
     * The index corresponding to the cell being compensated.
     */
    size_type index;

    /**
     * The integration limits for the cell being compensated.
     */
    DVector up;
    /**
     * The integration limits for the cell being compensated.
     */
    DVector lo;

  };

  /**
   * A vector (stack) of levels
   */
  typedef vector<Level> LevelVector;

  /**
   * The vector (stack) of levels
   */
  LevelVector levels;


  /**
   * This is a help struct to perform the divide-and-conquer slicing
   * of cells before starting the compensation procedure.
   */
  struct Slicer {

    /**
     * The constructor takes the number of dimensions of the function
     * approximated by the current cell, the ACDCGen object
     * controlling the generation and the lower-left and upper-right
     * corners of the cell to be sliced.
     */
    Slicer(DimType, ACDCGen &, const DVector &, const DVector &);

    /**
     * The constructor used internally when diagonally chopped-off
     * cells need to be sliced themselves.
     */
    Slicer(DimType Din, const Slicer & s, ACDCGenCell * cellin,
	   const DVector & loin, const DVector & xselin, const DVector & upin,
	   double fselin);

    /**
     * Destructor.
     */
    ~Slicer();

    /**
     * Called from both constructors to do the actual work.
     */
    void divideandconquer();

    /**
     * Initialize the procedure, finding the slicing points around the
     * current point
     */
    void init();

    /**
     * Do the slicing and increase the overestimate of the function in
     * the resulting cell. If a point with a higher function value has
     * been found repeat the slicing around that point etc.
     */
    void slice();

    /**
     * After slicing a cell, find the maximum function value found in
     * the resulting cell. Also set the minimum value found.
     */
    double shiftmaxmin();

    /**
     * Find the slice point of the current cell in the direction given.
     */
    void dohalf(DimType);

    /**
     * If split is in more than one dimensions check the overestimate
     * for the chopped-off cell.
     */
    void checkdiag(ACDCGenCell * cell, DimType d, double lod, double upd);

    /**
     * The dimension of the cell to be sliced.
     */
    DimType D;

    /**
     * The lower-left corner of the current cell.
     */
    DVector lo;
    /**
     * The  upper-right corner of the current cell.
     */
    DVector up;

    /**
     * The lower-left point found closest to the current
     * point which gives a function value below the overestimate.
     */
    DVector xcl;
    /**
     * The upper-right point found closest to the current point which
     * gives a function value below the overestimate.
     */
    DVector xcu;

    /**
     * The lower-left point furthest away from the
     * current point which gives a function value abov the
     * overestimate.
     */
    DVector xhl;
    /**
     * The upper-right point furthest away from the
     * current point which gives a function value abov the
     * overestimate.
     */
    DVector xhu;

    /**
     * The function values found for the xhl point.
     */
    DVector fhl;

    /**
     * The function values found for the xhu point.
     */
    DVector fhu;

    /**
     * The current point around which we are slicing.
     */
    DVector xsel;

    /**
     * The function value in the current point.
     */
    double fsel;

    /**
     * The current cell.
     */
    ACDCGenCell * current;

    /**
     * The cell which resulted from the first slicing procedure. This
     * is the first one to get an increased overestimate and is the
     * one to be compensated. All other cells with increased
     * overestimates are sub-cells to this one
     */
    ACDCGenCell * first;

    /**
     * The lower-left corner of the 'first' cell.
     */
    DVector firstlo;
    /**
     * The upper-right corner of the 'first' cell.
     */
    DVector firstup;

    /**
     * A pointer to the function to be used.
     */
    FncPtr f;

    /**
     * The epsilon() value obtained from the controlling
     * ACDCGen object.
     */
    double epsilon;

    /**
     * The margin() value obtained from the controlling
     * ACDCGen object.
     */
    double margin;

    /**
     * The dimensions to slice in rated by the resulting fractional
     * volume of the resulting slice. If the dimension is negative it
     * means that the cell should be slized from below.
     */
    multimap<double,DimType> rateslice;

    /**
     * The minimu function value found in the current sliced cell (set
     * by shiftmaxmin()).
     */
    double minf;

    /**
     * If true, then the whole original cell should compensated in the
     * continued generation.
     */
    bool wholecomp;

  };

public:

  /** The maximum recursion depth of the compensation so far. */
  static size_type maxsize;

private:

  /**
   * Copy constructor is private and not implemented.
   */
  ACDCGen(const ACDCGen &);

  /**
   * Assignment is private and not implemented.
   */
  ACDCGen & operator=(const ACDCGen &);

};

}

#include "ACDCGen.icc"

#endif
