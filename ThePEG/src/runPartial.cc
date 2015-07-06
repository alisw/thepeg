// -*- C++ -*-
//
// runPartial.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Utilities/DynamicLoader.h"
#include "ThePEG/Vectors/GenEventConverter.h"

int main(int argc, char * argv[]) {

  using namespace ThePEG;

  string run;

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    if ( arg == "-r" ) run = argv[++iarg];
    else if ( arg == "-l" ) DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    else if ( arg == "-L" ) DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    else if ( arg == "-h" ) {
    cerr << "Usage: " << argv[0]
	 << " [-l load-path] [-L first-load-path] run-file" << endl;
      return 3;
    }
    else
      run = arg;
  }

  if ( run.empty() ) {
    cerr << "No run-file specified." << endl;
    return 1;
  }

  try {

    PersistentIStream is(run);
    EGPtr eg;
    is >> eg;

    if ( !eg ) throw std::runtime_error("No generator found.");

    eg->initialize();
    // Open a persistent stream and read an OldEventGenerator object
    // from it and initialize the Event generator.


    for ( int i = 0; i < 10; ++i ) {

      EventPtr event = new_ptr(Event(PPair()));
      StepPtr firstStep = event->newStep();
      // Create an empty step.
      PPtr u = eg->getParticle(ParticleID::u);
      PPtr d = eg->getParticle(ParticleID::ubar);
      // Create a quark and an anti-quark.

      u->set3Momentum(Momentum3(ZERO, ZERO, 400.0*MeV));
      d->set3Momentum(Momentum3(ZERO, ZERO, -400.0*MeV));
      // Set the momentum of the quarks.
      u->antiColourNeighbour(d);
      d->colourNeighbour(u);
      // Setup their colour connections
      firstStep->addParticle(u);
      firstStep->addParticle(d);
      // Add the quarks to the first step

      PPtr p = eg->getParticle(ParticleID::pplus);
      p->set3Momentum(Momentum3(100.0*GeV, ZERO, ZERO));
      firstStep->addParticle(p);
      // Sometimes, the string fragmentation routine needs to shuffle some
      // energy around, so we add an extra proton for this purpose.

      event = eg->partialEvent(event);
      // Generate the event starting from the first step.

      //      event->removeParticle(p);
      // Remove the proton used for energy-momentum conservation.

      //      vector<PPtr> products;
      //      event->selectFinalState(inserter(products));
      // Extract the final state particles.

      cout << *event;
      // Do something interesting with the event.

      CLHEPMC::GenEvent * geneve = GenEventConverter::convert(*event);
      geneve->print(cout);
      delete geneve;

    }

    for ( int i = 0; i < 10; ++i ) {
      EventPtr event = eg->shoot();
      cout << *event;
      // Do something interesting with the event.

      using namespace CLHEPMC;

      CLHEPMC::GenEvent * geneve = GenEventConverter::convert(*event);
      geneve->print(cout);

      for ( GenEvent::particle_const_iterator it = geneve->particles_begin();
	    it != geneve->particles_end(); ++it ) {
	(**it).print();
      }

      delete geneve;
    }

    eg->finalize();
    // Tell the generator to write out statistics and stuff.

  }
  catch ( std::exception & e ) {
    cerr << e.what() << endl;
    return 1;
  }
  catch ( ... ) {
    cerr << "Unknown Exception\n";
    return 2;
  }

  return 0;
}

