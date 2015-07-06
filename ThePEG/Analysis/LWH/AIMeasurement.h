// -*- C++ -*-
//
// AIMeasurement.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIMeasurement_H
#define LWH_AIMeasurement_H

#ifndef LWH_USING_AIDA

/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IMeasurement {

public: 
    virtual ~IMeasurement() {}
    virtual double value() const = 0;
    virtual double errorPlus() const = 0;
    virtual double errorMinus() const = 0;
    virtual bool setValue(double value) = 0;
    virtual bool setErrorPlus(double errorPlus) = 0;
    virtual bool setErrorMinus(double errorMinus) = 0;
};

}

/** @endcond */

#else
#include "AIDA/IMeasurement.h"
#endif

#endif /* LWH_AIMeasurement_H */
