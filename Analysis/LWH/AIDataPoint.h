// -*- C++ -*-
//
// AIDataPoint.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIDataPoint_H
#define LWH_AIDataPoint_H



/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IMeasurement;

/**
 * Basic user-level interface class for holding and managing
 * a single set of "measurements".
 * 
 * @author The AIDA team (http://aida.freehep.org/)
 *
 */

class IDataPoint {

public: 
    virtual ~IDataPoint() {}
    virtual int dimension() const = 0;
    virtual IMeasurement * coordinate(int coord) = 0;
    virtual const IMeasurement * coordinate(int coord) const = 0;
};

}

/** @endcond */





#endif /* LWH_AIDataPoint_H */
