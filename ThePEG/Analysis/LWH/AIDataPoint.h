// -*- C++ -*-
//
// AIDataPoint.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIDataPoint_H
#define LWH_AIDataPoint_H

#ifndef LWH_USING_AIDA

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

#else
#include "AIDA/IDataPoint.h"
#endif

#endif /* LWH_AIDataPoint_H */
