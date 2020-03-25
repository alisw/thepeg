// -*- C++ -*-
//
// AIMeasurement.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIMeasurement_H
#define LWH_AIMeasurement_H



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





#endif /* LWH_AIMeasurement_H */
