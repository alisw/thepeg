// -*- C++ -*-
//
// Graphviz.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Graphviz class.
//

#include "GraphvizPlot.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

void GraphvizPlot::dofinish() {
  AnalysisHandler::dofinish();
  if (! _quiet )
    cout << "\nGraphvizPlot: plots can be generated like this:\n"
	 << "GraphvizPlot: 'dot -Tsvg " 
	 << generator()->filename() << '-'
	 << name() << '-'
	 << _eventNumber << ".dot > plot.svg'\n";
}

void GraphvizPlot::analyze(tEventPtr event, long, int, int) {
  if (event->number() != _eventNumber) return;

  // prepare dot file
  ostringstream fname;
  fname << generator()->filename() << '-' 
	<< name() << '-'
	<< event->number() << ".dot";
  ofstream hepmcdotfile(fname.str().c_str());

  printGraphviz(hepmcdotfile, event);

  hepmcdotfile.close();
}

void GraphvizPlot::persistentOutput(PersistentOStream & os) const {
  os << _eventNumber << _quiet;
}

void GraphvizPlot::persistentInput(PersistentIStream & is, int) {
  is >> _eventNumber >> _quiet;
}

ClassDescription<GraphvizPlot> GraphvizPlot::initGraphvizPlot;
// Definition of the static class description member.

void GraphvizPlot::Init() {

  static ClassDocumentation<GraphvizPlot> documentation
    ("There is no documentation for the GraphvizPlot class");

  static Parameter<GraphvizPlot,long> interfaceEventNumber
    ("EventNumber",
     "The number of the event that should be drawn.",
     &GraphvizPlot::_eventNumber, 1, 1, 1,
     false, false, Interface::lowerlim);
  interfaceEventNumber.setHasDefault(false);

  static Switch<GraphvizPlot,bool> interfaceQuiet
    ("Quiet",
     "Prevent GraphvizPlot from outputing instructions for how to generate "
     "the actual graph.",
     &GraphvizPlot::_quiet, false, true, false);
  static SwitchOption interfaceQuietVerbose
    (interfaceQuiet,
     "Verbose",
     "Allow output.",
     false);
  static SwitchOption interfaceQuietQuiet
    (interfaceQuiet,
     "Quiet",
     "Prevent output.",
     true);
  interfaceQuiet.setHasDefault(false);

}

