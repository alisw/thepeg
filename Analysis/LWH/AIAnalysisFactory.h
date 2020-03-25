// -*- C++ -*-
//
// AIAnalysisFactory.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIAnalysisFactory_H
#define LWH_AIAnalysisFactory_H



#include <string>

/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IDataPointSetFactory;
class IFitFactory;
class IFunctionFactory;
class IPlotterFactory;
class ITupleFactory;
class ITreeFactory;
class ITree;
class IHistogramFactory;

class IAnalysisFactory {

public:

  virtual ~IAnalysisFactory() {}

  virtual ITreeFactory * createTreeFactory() = 0;
  virtual IHistogramFactory * createHistogramFactory(ITree & tree) = 0;
  virtual IDataPointSetFactory * createDataPointSetFactory(ITree &) = 0;
  virtual ITupleFactory * createTupleFactory(ITree &) = 0;
  virtual IFunctionFactory * createFunctionFactory(ITree &) = 0;
  virtual IFitFactory * createFitFactory() = 0;
  virtual IPlotterFactory * createPlotterFactory(int = 0, char * * = 0,
						 const std::string & = "",
						 const std::string & = "") = 0;

};

}

/** @endcond */





#endif /* LWH_AIAnalysisFactory_H */
