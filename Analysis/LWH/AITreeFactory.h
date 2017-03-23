// -*- C++ -*-
//
// AITreeFactory.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AITreeFactory_H
#define LWH_AITreeFactory_H



#include <string>

/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class ITree;

class ITreeFactory {

public:

  virtual ~ITreeFactory() {}

  virtual ITree * create(const std::string &, const std::string & = "",
			 bool = false, bool = false,
			 const std::string & = "") = 0;

};

}

/** @endcond */





#endif /* LWH_AITreeFactory_H */
