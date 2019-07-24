// -*- C++ -*-
//
// ACDCGenCell.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ACDCGenCell_H
#define ACDCGenCell_H

#include "ACDCGenConfig.h"
#include "ACDCTraits.h"

namespace ACDCGenerator {

struct ACDCGenCellInfo;

/** ACDCGenCell is the class representing a generation cell in ACDCGen. */
class ACDCGenCell {

public:

  /**
   * Constructor taking the maximum value as argument.
   */
  inline ACDCGenCell(double newG);

  /**
   * Constructor taking the maximum value and the volume as argument.
   */
  inline ACDCGenCell(double newG, double newV);

  /**
   * The destuctor will also delete all child cells.
   */
  inline ~ACDCGenCell();

  /**
   * Choose a cell recursively according to their relative
   * overestimated integrals.
   * @param lo the lower-left corner of the chosen cell.
   * @param up the upper-right corner of the chosen cell.
   * @param rnd the random generator object used to throw dice to
   * choose sub-cell.
   * @return a pointer to the chosen cell.
   */
  template <typename RndType>
  inline ACDCGenCell * generate(DVector & lo, DVector & up, RndType * rnd);

  /**
   * Choose a cell recursively according to their relative
   * overestimated integrals.
   * @param lo the lower-left corner of the chosen cell.
   * @param up the upper-right corner of the chosen cell.
   * @param rndv a pre-generated set of random numbers (one for each
   * dimension) used to choose sub-cell and then rescales that random
   * number to be reused by the sub-cell.
   * @return a pointer to the chosen cell.
   */
  inline ACDCGenCell * generate(DVector & lo, DVector & up, DVector & rndv);

  /**
   * Find a cell. For a given phase space point, \a x, find the
   * corresponding cell. Afterwards, the \a up and \a lo vectors will
   * contain the upper-right and lower-left corners of the chosen
   * cell.
   */
  inline ACDCGenCell * getCell(DVector & lo, const DVector & x, DVector &  up);

  /**
   * Smooth out the levels. If one cell has an overestimated integral
   * which is less than \a frac of the adjacent one, rescale it to
   * avoid situations where it is never sampled.
   */
  inline void smooth(double frac);

  /**
   * Returns true if this cell has been split.
   */
  inline bool isSplit() const;

  /**
   * Recalculate (recursively) the overestimated integral for this
   * cell (and of the sub-cells) and return it. Optionally \a rescale
   * the overestimated integral.
   */
  inline double doMaxInt(double rescale = 1.0);

  /**
   * Return the last calculated the overestimated integral for this
   * cell.
   */
  inline double maxInt() const;

  /**
   * Split this cell into two. The cell is split along the \a newDim
   * direction, where the lower and upper limit is given by \a lo and
   * \a up and the point of division is given by \a newDiv.
   */
  inline void splitme(double lo, double newDiv, double up, DimType newDim);

  /**
   * Set a new overestimated maximum function value in this cell.
   */
  inline void g(double newG);

  /**
   * Return the overestimated maximum function value in this cell.
   */
  inline double g() const;

  /**
   * Return the volume of this cell.
   */
  inline double v() const;

  /**
   * Return the direction in which it has been split. Return -1 if it
   * has not been split.
   */
  inline DimType dim() const;

  /**
   * Return the point of division in the dim() direction. Return -1.0
   * if it has not been split.
   */
  inline double div() const;

  /**
   * Return the upper sub-cell. Return null if it has not been split.
   */
  inline ACDCGenCell * upper() const;

  /**
   * Return the lower sub-cell. Return null if it has not been split.
   */
  inline ACDCGenCell * lower() const;

  /**
   * Return the number of cells in this sub-tree which have not been
   * split.
   */
  inline int nBins() const;

  /**
   * Return the maximum depth of this sub-tree.
   */
  inline int depth() const;

  /**
   * Append ACDCGenCellInfo objects describing this subtree to a given
   * vector.
   * @param lo the lower-left corner of this cell.
   * @param up the upper-right corner of this cell.
   * @param v the vector where the ACDCGenCellInfo objects will be appended.
   */
  inline void extract(DVector & lo, DVector & up,
		      vector<ACDCGenCellInfo> & v) const;

  /**
   * Get the index of the given cell.
   */
  inline long getIndex(const ACDCGenCell * c) const;

  /**
   * Helper function for getIndex(const ACDCGenCell *) with an extra
   * argument to use as counter.
   */
  inline long getIndex(const ACDCGenCell * c, long & indx) const;
   
  /**
   * Return the cell corresponding to the given index \a i.
   */
  inline ACDCGenCell * getCell(long i);

  /**
   * Helper function for getCell(long) with an extra argument to use as
   * counter.
   */
  inline ACDCGenCell * getCell(long i, long & indx);

public:

  /**
   * If the cell has not been split this is the overestimated maximum
   * function value in this cell. Otherwise it is the weighted
   * average of the sub-cells values.
   */
  double theG;

  /**
   * The volume of this cell.
   */
  double theV;

  /**
   * Pointers to the upper sub-cell.
   */
  ACDCGenCell * theUpper;

  /**
   * Pointers to the lower sub-cell.
   */
  ACDCGenCell * theLower;

  /**
   * The point of division in the dim() direction.
   */
  double theDivision;

  /**
   * The direction in which it has been split.
   */
  DimType theSplitDimension;

private:

  /**
   * Default constructor is private and not implemented.
   */
  ACDCGenCell();

  /**
   * Copy constructor is private and not implemented.
   */
  ACDCGenCell(const ACDCGenCell &);

  /**
   * Assignment is private and not implemented.
   */
  ACDCGenCell & operator=(const ACDCGenCell &) = delete;

};


/**
 * This is a class describing cells to the outside world to be used
 * for debugging purposes. They only make sense if extracted with the
 * ACDCGenCell::extract function.
 */
struct ACDCGenCellInfo {

  /** the integer used for indices. */
  typedef vector<ACDCGenCellInfo>::size_type Index;

  /**
   * The overestimated maximum function value of this cell.
   */
  double g;

  /**
   * The volume of the corresponding cell.
   */
  double v;

  /**
   * The upper-right corner of the corresponding cell.
   */
  DVector up;
  /**
   * The lower-left corner of the corresponding cell.
   */
  DVector lo;

  /**
   * The index of the upper sub-cells in the vector in which the
   * corresponding cell was inserted.
   */
  Index iup;

  /**
   * The index of the lower sub-cell in the vector in which the
   * corresponding cell was inserted.
   */
  Index ilo;

};

}

#include "ACDCGenCell.icc"

#endif
