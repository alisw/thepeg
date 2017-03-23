// -*- C++ -*-
//
// AIDataPointSetFactory.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIDataPointSetFactory_H
#define LWH_AIDataPointSetFactory_H



/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IDataPointSet;
class IHistogram1D;

class IDataPointSetFactory {

public:
    virtual ~IDataPointSetFactory() { /* nop */; }
    virtual IDataPointSet *
    create(const std::string &, const std::string &, int ) = 0;
    virtual IDataPointSet * create(const std::string &, int) = 0;
    virtual IDataPointSet *
    createY(const std::string &, const std::string &,
	    const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createY(const std::string &, const std::string &,
	    const std::vector<double> &, const std::vector<double> &,
	    const std::vector<double>  &) = 0;
    virtual IDataPointSet *
    createY(const std::string &, const std::vector<double> &,
	    const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createY(const std::string &, const std::vector<double> &,
	    const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createX(const std::string &, const std::string &,
	    const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createX(const std::string &, const std::string &,
	    const std::vector<double> &, const std::vector<double> &,
	    const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createX(const std::string &, const std::vector<double> &,
	    const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createX(const std::string &, const std::vector<double> &,
	    const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXY(const std::string &, const std::string &,
	     const std::vector<double> &, const std::vector<double> &,
	     const std::vector<double> &, const std::vector<double> &,
	     const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXY(const std::string &, const std::string &,
	     const std::vector<double> &, const std::vector<double> &,
	     const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXY(const std::string &, const std::vector<double> &,
	     const std::vector<double> &, const std::vector<double> &,
	     const std::vector<double> &, const std::vector<double> &,
	     const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXY(const std::string &, const std::vector<double> &,
	     const std::vector<double> &, const std::vector<double> &,
	     const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXYZ(const std::string &, const std::string &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXYZ(const std::string &, const std::string &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXYZ(const std::string &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createXYZ(const std::string &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &, const std::vector<double> &,
	      const std::vector<double> &) = 0;
    virtual IDataPointSet *
    createCopy(const std::string &, const IDataPointSet &) = 0;
    virtual bool destroy(IDataPointSet *) = 0;
    virtual IDataPointSet * create(const std::string &, const IHistogram1D &,
				   const std::string & = "") = 0;

};

}

/** @endcond */





#endif /* LWH_AIDataPointSetFactory_H */
