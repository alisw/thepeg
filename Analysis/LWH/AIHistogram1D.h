// -*- C++ -*-
//
// AIHistogram1D.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIHistorgram1D_H
#define LWH_AIHistorgram1D_H



/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IAnnotation;

class IBaseHistogram {

public: 

  virtual ~IBaseHistogram() {}

  virtual std::string title() const = 0;
  virtual bool setTitle(const std::string & title) = 0;
  virtual int dimension() const = 0;
  virtual bool reset() = 0;
  virtual int entries() const = 0;

};

class IHistogram : virtual public IBaseHistogram {

public: 

  virtual ~IHistogram() {}

  virtual int allEntries() const = 0;
  virtual int extraEntries() const = 0;
  virtual double equivalentBinEntries() const = 0;
  virtual double sumBinHeights() const = 0;
  virtual double sumAllBinHeights() const = 0;
  virtual double sumExtraBinHeights() const = 0;
  virtual double minBinHeight() const = 0;
  virtual double maxBinHeight() const = 0;

};

class IAxis;

class IHistogram1D: virtual public IHistogram {

public:

  virtual ~IHistogram1D() {}

  virtual bool fill(double x, double weight = 1.) = 0;
  virtual double binMean(int index) const = 0;
  virtual int binEntries(int index) const = 0;
  virtual double binHeight(int index) const = 0;
  virtual double binError(int index) const = 0;
  virtual double mean() const = 0;
  virtual double rms() const = 0;
  virtual const IAxis & axis() const = 0;
  virtual int coordToIndex(double coord) const = 0;
  virtual bool add(const IHistogram1D & hist) = 0;
  virtual bool scale(double scaleFactor) = 0;

};

}

/** @endcond */





#endif /* LWH_AIHistorgram1D_H */
