// -*- C++ -*-
//
// runThePEG.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/DynamicLoader.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Repository/Main.h"
#include "ThePEG/Repository/Repository.h"
#include <config.h>

int main(int argc, char * argv[]) {
  using namespace ThePEG;

  string run;
  long N = -1;
  long seed = 0;
  string mainclass;
  bool tics = false;
  bool resume = false;
  string tag = "";
  string setupfile = "";

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    if ( arg == "-r" ) run = argv[++iarg];
    else if ( arg == "-x" ) mainclass = argv[++iarg];
    else if ( arg == "-m" ) setupfile = argv[++iarg];
    else if ( arg == "-s" ) DynamicLoader::load(argv[++iarg]);
    else if ( arg.substr(0,2) == "-s" )
      DynamicLoader::load(arg.substr(2));
    else if ( arg == "-l" ) DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    else if ( arg == "-L" ) DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    else if ( arg == "-d" ) Debug::setDebug(atoi(argv[++iarg]));
    else if ( arg.substr(0,2) == "-d" )
      Debug::setDebug(atoi(arg.substr(2).c_str()));
    else if ( arg.substr(0,2) == "-D" ) DebugItem::setDebugItem(arg.substr(2));
    else if ( arg == "-N" ) N = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-N" ) N = atoi(arg.substr(2).c_str());
    else if ( arg == "--seed" || arg == "-seed" ) seed = atol(argv[++iarg]);
    else if ( arg == "--tics" || arg == "-tics" ) tics = true;
    else if ( arg == "--resume" ) resume = true;
    else if ( arg == "-t" ) tag = argv[++iarg];
    else if ( arg.substr(0,2) == "-t" ) tag = arg.substr(2);
    else if ( arg.substr(0,6) == "--tag=" ) tag = arg.substr(6);
    else if ( arg == "--help" || arg == "-h" ) {
    cerr << "Usage: " << argv[0] << " [-d {debuglevel|-debugitem}] "
	 << "[-l load-path] [-L first-load-path] [-m setup-file] run-file" << endl;
      return 3;
    }
    else if ( arg == "-v" || arg == "--version" ) {
      cout << PACKAGE_VERSION << endl;
      return 0;
    }
    else
      run = arg;
  }

  if ( Debug::level ) Debug::unmaskFpuErrors();

  if ( run.empty() ) {
    cerr << "No run-file specified." << endl;
    return 1;
  }

  try {

    EGPtr eg;
    if ( run == "-" ) {
      PersistentIStream is(cin);
      is >> eg;
    } else {
      PersistentIStream is(run);
      is >> eg;
    }

    breakThePEG();

    if ( !eg ) {
      cout << "Could not find or read the requested EventGenerator." << endl;
      return 1;
    }

    if ( setupfile.size() ) {
      string msg = Repository::modifyEventGenerator(*eg, setupfile, cout);
      if ( ! msg.empty() ) cerr << msg << '\n';
    }

    if ( seed > 0 ) eg->setSeed(seed);
    if ( !tag.empty() ) eg->addTag(tag);
    if ( !mainclass.empty() ) {
      Main::arguments(vector<string>(argv + 1, argv + argc));
      Main::N(N);
      if ( !eg->loadMain(mainclass) )
	std::cout << "Main class file '" << mainclass << "' not found." << endl;
    } else {
      eg->go(resume? -1: 1, N, tics);
    }
  }
  catch ( Exception & e ) {
    cerr << "Unexpected exception caught: " << e.what() << endl;
    e.handle();
    return 1;
  }
  catch ( std::exception & e ) {
    cerr << "Unexpected exception caught: " << e.what() << endl;
    return 1;
  }
  catch ( ... ) {
    breakThePEG();
    cerr << "Unknown Exception\n";
    return 2;
  }

  return 0;
}

