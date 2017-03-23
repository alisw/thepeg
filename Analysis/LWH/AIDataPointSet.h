// -*- C++ -*-
//
// AIDataPointSet.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIDataPointSet_H
#define LWH_AIDataPointSet_H



/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IAnnotation;
class IDataPoint;

/**
 * Basic user-level interface class for holding and managing
 * a single set of "data points".
 *
 * @author The AIDA team (http://aida.freehep.org/)
 *
 */

class IDataPointSet {

public: 
    virtual ~IDataPointSet() { /* nop */; }
    virtual IAnnotation & annotation() = 0;
    virtual const IAnnotation & annotation() const = 0;
    virtual std::string title() const = 0;
    virtual bool setTitle(const std::string & title) = 0;
    virtual int dimension() const = 0;
    virtual void clear() = 0;
    virtual int size() const = 0;
    virtual IDataPoint * point(int index) = 0;
    virtual bool setCoordinate(int coord, const std::vector<double>  & val, const std::vector<double>  & err) = 0;
    virtual bool setCoordinate(int coord, const std::vector<double>  & val, const std::vector<double>  & errp, const std::vector<double>  & errm) = 0;
    virtual const IDataPoint * point(int index) const = 0;
    virtual IDataPoint * addPoint() = 0;
    virtual bool addPoint(const IDataPoint & point) = 0;
    virtual bool removePoint(int index) = 0;
    virtual double lowerExtent(int coord) const = 0;
    virtual double upperExtent(int coord) const = 0;
    virtual bool scale(double scaleFactor) = 0;
    virtual bool scaleValues(double scaleFactor) = 0;
    virtual bool scaleErrors(double scaleFactor) = 0;
    virtual void * cast(const std::string & className) const = 0;
};

}

/** @endcond */





#endif /* LWH_AIDataPointSet_H */
