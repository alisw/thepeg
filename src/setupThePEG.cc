// -*- C++ -*-
//
// setupThePEG.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

// macro is passed in from -D compile flag
#ifndef THEPEG_PKGLIBDIR
#error Makefile.am needs to define THEPEG_PKGLIBDIR
#endif

#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Utilities/DynamicLoader.h"

int main(int argc, char * argv[]) {
  using namespace ThePEG;

  Debug::level = 1;

  // macro is passed in from -D compile flag
  string repo = string(THEPEG_PKGLIBDIR) + "/ThePEGDefaults.rpo";

  string repout;
  string file;
  bool init = false;
  vector<string> globlib;
  vector<string> preread;
  vector<string> appread;

  Repository repository;

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    if ( arg == "-d" ) Debug::setDebug(atoi(argv[++iarg]));
    else if ( arg.substr(0,2) == "-d" )
      Debug::setDebug(atoi(arg.substr(2).c_str()));
    else if ( arg == "-r" ) repo = argv[++iarg];
    else if ( arg == "-o" ) repout = argv[++iarg];
    else if ( arg == "--init" || arg == "-init" ) {
      init = true;
      Debug::level = 0;
    }
    else if ( arg == "--exitonerror" ) repository.exitOnError() = 1;
    else if ( arg == "-s" ) {
      DynamicLoader::load(argv[++iarg]);
      repository.globalLibraries().push_back(argv[iarg]);
      globlib.push_back(argv[iarg]);
    }
    else if ( arg.substr(0,2) == "-s" ) {
      DynamicLoader::load(arg.substr(2));
      repository.globalLibraries().push_back(arg.substr(2));
      globlib.push_back(arg.substr(2));
    }
    else if ( arg == "-l" ) DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    else if ( arg == "-L" ) DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    else if ( arg == "-i" ) {
      Repository::appendReadDir(argv[++iarg]);
      appread.push_back(argv[iarg]);
    }
    else if ( arg.substr(0,2) == "-i" ) {
      Repository::appendReadDir(arg.substr(2));
      appread.push_back(arg.substr(2));
    }
    else if ( arg == "-I" ) {
      Repository::prependReadDir(argv[++iarg]);
      preread.push_back(argv[iarg]);
    }
    else if ( arg.substr(0,2) == "-I" ) {
      Repository::prependReadDir(arg.substr(2));
      preread.push_back(arg.substr(2));
    }
    else if ( arg == "-h" || arg == "--help" ) {
      cerr << "Usage: " << argv[0]
	 << " {cmdfile} [-d {debuglevel|-debugitem}] [-r input-repository-file]"
	 << " [-l load-path] [-L first-load-path]" << endl;
      return 3;
    }
    else if ( arg == "-v" || arg == "--version" ) {
      cout << Repository::version() << endl;
      return 0;
    }
    else
      file = arg;
  }

  if ( Debug::level ) Debug::unmaskFpuErrors();

  try {

    if ( init ) {
      breakThePEG();
      if ( repout.empty() ) repout = repo;
      else {
	string msg = repository.load(repo);
	if ( ! msg.empty() ) cerr << msg << '\n';
	for ( unsigned int i = 0; i < globlib.size(); ++i )
	  repository.globalLibraries().push_back(globlib[i]);
	for ( unsigned int i = 0; i < appread.size(); ++i )
	  Repository::appendReadDir(appread[i]);
	for ( unsigned int i = 0; i < preread.size(); ++i )
	  Repository::prependReadDir(preread[i]);
      }
      {
	HoldFlag<> setup(InterfaceBase::NoReadOnly);
	if ( file.empty() ) file = "ThePEGDefaults.in";
	string msg = repository.read(file, cout);
	if ( ! msg.empty() ) cerr << msg << '\n';
	repository.update();
      }
      repository.save(repout);
    } else {
      string msg = repository.load(repo);
      if ( ! msg.empty() ) cerr << msg << '\n';
      for ( unsigned int i = 0; i < globlib.size(); ++i )
	repository.globalLibraries().push_back(globlib[i]);
      for ( unsigned int i = 0; i < appread.size(); ++i )
	Repository::appendReadDir(appread[i]);
      for ( unsigned int i = 0; i < preread.size(); ++i )
	Repository::prependReadDir(preread[i]);
      breakThePEG();
      if ( file.size() && file != "-" ) {
	if ( file == "--java" || file == "-java" )
	  repository.read(cin, cout, "-*-ready-*-\n");
	else {
	  string msg = repository.read(file, cout);
	  if ( ! msg.empty() ) cerr << msg << '\n';
	}
      } else {
	repository.read(cin, cout, "ThePEG> ");
      }
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
    cerr << "Unknown Exception\n";
    return 2;
  }

  return 0;
}

