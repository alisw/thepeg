// -*- C++ -*-
//
// runEventLoop.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DynamicLoader.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Vectors/HepMCConverter.h"
#include "HepMC/GenEvent.h"

namespace ThePEG {

template<>
struct HepMCTraits<HepMC::GenEvent>:
    public HepMCTraitsBase<HepMC::GenEvent,HepMC::GenParticle,
                           HepMC::GenVertex,HepMC::Polarization,
			   HepMC::PdfInfo> {};

}

int main(int argc, char * argv[]) {
  using namespace ThePEG;

  // First we get some command-line arguments

  string run;
  long N = -1;
  long seed = 0;

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    // Specifying a run file
    if ( arg == "-r" ) run = argv[++iarg];
    // Append a path for the dynamic loader
    else if ( arg == "-l" ) DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    // Prepend a path for the dynamic loader
    else if ( arg == "-L" ) DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    // Set debug level
    else if ( arg == "-d" ) Debug::setDebug(atoi(argv[++iarg]));
    else if ( arg.substr(0,2) == "-d" )
      Debug::setDebug(atoi(arg.substr(2).c_str()));
    // Set number of events
    else if ( arg == "-N" ) N = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-N" ) N = atoi(arg.substr(2).c_str());
    // Set random seed
    else if ( arg == "-seed" ) seed = atoi(argv[++iarg]);
    // Print (out of date) help message
    else if ( arg == "-h" ) {
    cerr << "Usage: " << argv[0] << " [-d {debuglevel|-debugitem}] "
	 << "[-l load-path] [-L first-load-path] run-file" << endl;
      return 3;
    }
    else
      // Any other argument is treated as a run file
      run = arg;
  }

  if ( Debug::level ) Debug::unmaskFpuErrors();

  if ( run.empty() ) {
    cerr << "No run-file specified." << endl;
    return 1;
  }

  try {

    // Create a persistent stream and read in an OldEventGenerator from
    // the run file
    PersistentIStream is(run);
    EGPtr eg;
    is >> eg;

    breakThePEG();

    if ( eg ) {

      if ( seed > 0 ) eg->setSeed(seed);

      // Initialize the Event generator
      eg->initialize();

      // Get number of events
      if ( N < 0 ) N = eg->N();

      // HERE IS THE MAIN EVENT LOOP
      for ( int ieve = 0; ieve < N; ++ieve ) {

	// Generate an event 
	EventPtr event = eg->shoot();

	// Convert to a HepMC::GenEvent
	HepMC::GenEvent * geneve =
	  HepMCConverter<HepMC::GenEvent>::convert(*event);

	// Do whatever you want with the event here
	if ( ieve < 10 ) geneve->print(cout);

	// Don't forget to delete the HepMC::GenEvent (The
	// ThePEG::Event is automatically garbage collected)
	delete geneve;

      }

      // End the run
      eg->finalize();

    } else std::cout << "eg = nil" << endl;

  }
  catch ( std::exception & e ) {
    cerr << e.what() << endl;
    return 1;
  }
  catch ( ... ) {
    breakThePEG();
    cerr << "Unknown Exception\n";
    return 2;
  }

  return 0;
}

