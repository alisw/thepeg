// -*- C++ -*-
//
// MultiEventGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MultiEventGenerator class.
//

#include "MultiEventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include <ctime>

using namespace ThePEG;

MultiEventGenerator::~MultiEventGenerator() {}

IBPtr MultiEventGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr MultiEventGenerator::fullclone() const {
  return new_ptr(*this);
}

string MultiEventGenerator::removeInterface(string cmd) {
  string noun = StringUtils::car(cmd);
  IBPtr ip = BaseRepository::getObjectFromNoun(noun);
  const InterfaceBase * ifb = BaseRepository::
    FindInterface(ip, BaseRepository::getInterfaceFromNoun(noun));
  string posarg = BaseRepository::getPosArgFromNoun(noun);

  for ( string::size_type i = 0; i < theObjects.size(); ++i ) {
    if ( theObjects[i] == ip && theInterfaces[i] == ifb->name() &&
	 thePosArgs[i] == posarg ) {
      theObjects.erase(theObjects.begin() + i);
      theInterfaces.erase(theInterfaces.begin() + i);
      thePosArgs.erase(thePosArgs.begin() + i);
      theValues.erase(theValues.begin() + i);
      return "";
    }
  }
  return "No such object/interface defined for this MultiEventGenerator.";
}

string MultiEventGenerator::addInterface(string cmd) {
  return addInterface(cmd, false);
}

string MultiEventGenerator::addRndInterface(string cmd) {
  return addInterface(cmd, true);
}

string MultiEventGenerator::addInterface(string cmd, bool rnd) {
  breakThePEG();
  string noun = StringUtils::car(cmd);
  IBPtr ip = BaseRepository::getObjectFromNoun(noun);
  const InterfaceBase * ifb = BaseRepository::
    FindInterface(ip, BaseRepository::getInterfaceFromNoun(noun));
  string posarg = BaseRepository::getPosArgFromNoun(noun);
  cmd = StringUtils::cdr(cmd);
  if ( cmd.empty() ) return "Error: empty argument list.";

  string ret;
  string oldvalue = ifb->exec(*ip, "get", posarg);
  StringVector args;
  try {
    if ( rnd ) {
      do {
	args.push_back(StringUtils::car(cmd));
	cmd = StringUtils::cdr(cmd);
      } while ( !cmd.empty() );
      if ( args.size() < 3 ) return "Error: Argument list should be 'N min max mean width'.";
      int N = atoi(args[0].c_str());
      string vmin = args[1];
      string vmax = args[2];
      string vmean = "0";
      if ( args.size() > 3 ) vmean = args[3];
      string vwidth = "0";
      if ( args.size() > 4 ) vwidth = args[4];
      string arg = "RND " + vmin + " " + vmax + " " + vmean + " " + vwidth;
      args = vector<string>(N, arg);
      ifb->exec(*ip, "set", vmin);
      ifb->exec(*ip, "set", vmax);
      ifb->exec(*ip, "set", posarg + " " + oldvalue);
    } else {
      do {
	args.push_back(StringUtils::car(cmd, ","));
	cmd = StringUtils::cdr(cmd, ",");
      } while ( !cmd.empty() );
      for ( string::size_type i = 0; i < args.size(); ++i )
	ifb->exec(*ip, "set", args[i]);
    }
  }
  catch (const Exception & e) {
    e.handle();
    ret = "Error: " + e.message();
  }
  ifb->exec(*ip, "set", posarg + " " + oldvalue);
  if ( !ret.empty() ) return ret;

  for ( string::size_type i = 0; i < theObjects.size(); ++i ) {
    if ( theObjects[i] == ip && theInterfaces[i] == ifb->name() &&
	 thePosArgs[i] == posarg ) {
      if ( rnd || theValues[i][0].substr(0,3) == "RND" ) theValues[i] = args;
      else theValues[i].insert(theValues[i].end(), args.begin(), args.end());
      return "";
    }
  }

  theObjects.push_back(ip);
  theInterfaces.push_back(ifb->name());
  thePosArgs.push_back(posarg);
  theValues.push_back(args);
  return "";
}

void MultiEventGenerator::addTag(string tag) {
  if ( tag[0] == '#' ) {
    string::size_type dash = tag.find('-');
    if ( dash == string::npos )
      firstSubrun = lastSubrun = atoi(tag.substr(1).c_str());
    else {
      firstSubrun = atoi(tag.substr(1, dash - 1).c_str());
      lastSubrun = atoi(tag.substr(dash + 1).c_str());
    }
  }
  EventGenerator::addTag(tag);
}



void MultiEventGenerator::doGo(long next, long maxevent, bool tics) {

  if ( theObjects.empty() || next < 0 ) {
    EventGenerator::doGo(next, maxevent, tics);
    return;
  }

  if ( maxevent >= 0 ) N(maxevent);

  vector<const InterfaceBase *> interfaces;
  long nargs = 1;
  for ( string::size_type i = 0; i < theObjects.size(); ++i ) {
    nargs *= theValues[i].size();
    interfaces.push_back(BaseRepository::FindInterface(theObjects[i],
						       theInterfaces[i]));
  }
  if ( theSeparateRandom ) {
    theSeparateRandom->init();
    theSeparateRandom->initrun();
    const InterfaceBase * ifb = BaseRepository::FindInterface(theSeparateRandom, "Seed");
    ifb->exec(*theSeparateRandom, "set", "0"); 
  }

  openOutputFiles();

  string baseName = runName();

  if ( tics ) tic(next - 1, nargs*N());
  for ( long iargs = 0; iargs < nargs; ++iargs ) {

    ostringstream subname;
    subname << baseName << ":" << iargs + 1;
    runName(subname.str());

    string head = heading(iargs, interfaces, baseName);

    if ( ( firstSubrun > 0 && iargs + 1 < firstSubrun ) ||
	 ( lastSubrun > 0 && iargs + 1 > lastSubrun ) ) {
      if ( theSeparateRandom ) {
	// This is needed to ensure the same random settings for a
	// given sub-run irrespectively if previous sub-runs have been
	// included or not.
	theSeparateRandom->reset();
	theSeparateRandom->init();
	theSeparateRandom->initrun();
      }
      continue;
    }
      

    log() << head;
    out() << head;

    reset();
    for_each(objects(), mem_fun(&InterfacedBase::reset));
    
    init();
    initrun();

    ieve = next-1;

    try {
      while ( shoot() ) {
	if ( tics ) tic(ieve + iargs*N(), nargs*N());
      }
    }
    catch ( ... ) {
      finish();
      throw;
    }
    finish();

  }

  runName(baseName);

  finally();

}

string MultiEventGenerator::
heading(long iargs, const vector<const InterfaceBase *> & interfaces,
	string baseName) const {
  ostringstream os;
  long div = 1;
  if ( iargs > 0 ) os << endl;
      
  os << ">> " << baseName << " sub-run number " << iargs + 1
     << " using the following interface values:" << endl;

  for ( string::size_type i = 0; i < theObjects.size(); ++i ) {
    long iarg = (iargs/div)%theValues[i].size();
    string sval = theValues[i][iarg];
    if ( theValues[i][iarg].substr(0,3) == "RND" ) {
      double vmin, vmax, vmean, vwidth;
      istringstream is(theValues[i][iarg].substr(3));
      is >> vmin >> vmax >> vmean >> vwidth;
      double val = randomArg().rnd(vmin, vmax);
      if ( vwidth > 0.0 ) do {
	  val = randomArg().rndGauss(vwidth, vmean);
	} while ( val < vmin || val > vmax );
      ostringstream ssv;
      ssv << val;
      sval = ssv.str();
    }
    interfaces[i]->exec(*theObjects[i], "set",
			thePosArgs[i] + " " + sval);
    os << "   set " << theObjects[i]->name() << ":" << theInterfaces[i];
    if ( !thePosArgs[i].empty() ) os << "[" << thePosArgs[i] << "]";
    os << " " << sval << endl;
    div *= theValues[i].size();
  }
  os << endl;
  return os.str();
}

void MultiEventGenerator::persistentOutput(PersistentOStream & os) const {
  os << theObjects << theInterfaces << thePosArgs << theValues
     << firstSubrun << lastSubrun << theSeparateRandom;
}

void MultiEventGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theObjects >> theInterfaces >> thePosArgs >> theValues
     >> firstSubrun >> lastSubrun >> theSeparateRandom;
}

IVector MultiEventGenerator::getReferences() {
  IVector ret = EventGenerator::getReferences();
  ret.insert(ret.end(), theObjects.begin(), theObjects.end());
  return ret;
}

void MultiEventGenerator::rebind(const TranslationMap & trans)
  {
  for ( string::size_type i = 0; i < theObjects.size(); ++i )
    theObjects[i] = trans.translate(theObjects[i]);
  EventGenerator::rebind(trans);
}

ClassDescription<MultiEventGenerator>
MultiEventGenerator::initMultiEventGenerator;
// Definition of the static class description member.

void MultiEventGenerator::Init() {

  static ClassDocumentation<MultiEventGenerator> documentation
    ("The ThePEG::MultiEventGenerator class is derived from the "
     "ThePEG::EventGenerator and is capable of making "
     "several runs with a pre-defined set of parameter and switch values.");

  static Command<MultiEventGenerator> interfaceAddInterface
    ("AddInterface",
     "If arguments are given on the form 'object-name:interface-name arg1, "
     "arg2, arg3' or 'object-name:vectorinterface-name[pos] arg1, arg2, arg3' "
     "the generator will be run three times with the corresonding interface of "
     "the given object set to arg1, arg2, arg3 in each run respectively. If "
     "another interface with e.g. 4 different arguments, the generator will "
     "be run 12 times once for each combination of arguments. If called with "
     "an object and interface wich has already been given in a previous call, "
     "the new arguments will be added to the previously specified list without "
     "checking if any argument is doubled.",
     &MultiEventGenerator::addInterface);

  static Command<MultiEventGenerator> interfaceAddRndInterface
    ("AddRndInterface",
     "If arguments are given on the form 'object-name:interface-name N min max "
     "mean width'"" or 'object-name:vectorinterface-name[pos] N min max mean width' "
     "the generator will be run N times with the corresonding interface of "
     "the given object set to a random value between min and max according to a "
     "Gaussian distribution with the given mean and width (if the width is absent "
     "or zero a flat distribution between min and max will be used instead) . If "
     "another interface with e.g. 4 different arguments, the generator will "
     "be run N*4 times once for each combination of arguments and the specified "
     "interface will get a new random value each time.. If called with "
     "an object and interface wich has already been given in a previous call to "
     "<interface>AddInterface</interface> or <interface>AddRndInterface</interface> "
     "the previous call will be ignored.",
     &MultiEventGenerator::addRndInterface);

  static Command<MultiEventGenerator> interfaceRemoveInterface
    ("RemoveInterface",
     "If arguments are given on the form 'object-name:interface-name' and "
     "the same interface and object was previously with an "
     "<interface>AddInterface</interface>}, the corresponding arguments are "
     "removed and the interfaced will be left unchanged during the generation.",
     &MultiEventGenerator::removeInterface);


  static Reference<MultiEventGenerator,RandomGenerator> interfaceSeparateRandom
    ("SeparateRandom",
     "A separate random number generator used for <interface>AddRndInterface</interface>"
     " to ensure reproducible sequences of interface values. If null, the standard "
     "random generator will be used instead.",
     &MultiEventGenerator::theSeparateRandom, true, false, true, true, false);

  interfaceAddInterface.rank(10.7);
  interfaceRemoveInterface.rank(10.5);

}

