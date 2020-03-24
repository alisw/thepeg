// -*- C++ -*-
//
// Repository.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Repository class.
//

// macro is passed in from -D compile flag
#ifndef THEPEG_PKGLIBDIR
#error Makefile.am needs to define THEPEG_PKGLIBDIR
#endif

#include "Repository.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Repository/Strategy.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Utilities/DynamicLoader.h"
#include "ThePEG/Utilities/StringUtils.h"

#include <iterator>
#include <chrono>

#include <config.h>

// readline options taken from
// http://autoconf-archive.cryp.to/vl_lib_readline.html 
// Copyright © 2008 Ville Laurikari <vl@iki.fi> 
// Copying and distribution of this file, with or without
// modification, are permitted in any medium without royalty provided
// the copyright notice and this notice are preserved.

#ifdef HAVE_LIBREADLINE
#  if defined(HAVE_READLINE_READLINE_H)
#    include <readline/readline.h>
#  elif defined(HAVE_READLINE_H)
#    include <readline.h>
#  else
extern "C" char *readline (const char *);
#  endif
#endif

#ifdef HAVE_READLINE_HISTORY
#  if defined(HAVE_READLINE_HISTORY_H)
#    include <readline/history.h>
#  elif defined(HAVE_HISTORY_H)
#    include <history.h>
#  else
extern "C" void add_history (const char *);
#  endif
#endif

using namespace ThePEG;

ParticleMap & Repository::defaultParticles() {
  static ParticleMap theMap;
  return theMap;
}


ParticleDataSet & Repository::particles() {
  static ParticleDataSet theSet;
  return theSet;
}

MatcherSet & Repository::matchers() {
  static MatcherSet theSet;
  return theSet;
}

Repository::GeneratorMap & Repository::generators() {
  static GeneratorMap theMap;;
  return theMap;
}

string & Repository::currentFileName() {
  static string theCurrentFileName;
  return theCurrentFileName;
}

int & Repository::exitOnError() {
  static int exitonerror = 0;
  return exitonerror;
}

void Repository::cleanup() {
  generators().clear();
}

void Repository::Register(IBPtr ip) {
  BaseRepository::Register(ip);
  registerParticle(dynamic_ptr_cast<PDPtr>(ip));
  registerMatcher(dynamic_ptr_cast<PMPtr>(ip));
}

void Repository::Register(IBPtr ip, string newName) {
  DirectoryAppend(newName);
  BaseRepository::Register(ip, newName);
  registerParticle(dynamic_ptr_cast<PDPtr>(ip));
  registerMatcher(dynamic_ptr_cast<PMPtr>(ip));
}

void Repository::registerParticle(tPDPtr pd) {
  if ( !pd ) return;
  if ( !member(particles(), pd) ) {
    particles().insert(pd);
    CreateDirectory(pd->fullName());
  }
  if ( pd->id() == 0 ) return;
  if ( !member(defaultParticles(), pd->id()) )
    defaultParticles()[pd->id()] = pd;
  for ( MatcherSet::iterator it = matchers().begin();
	it != matchers().end(); ++it) (*it)->addPIfMatch(pd);
}

void Repository::registerMatcher(tPMPtr pm) {
  if ( !pm || member(matchers(), pm) ) return;
  pm->addPIfMatchFrom(particles());
  for ( MatcherSet::iterator it = matchers().begin();
	it != matchers().end(); ++it) {
    (*it)->addMIfMatch(pm);
    pm->addMIfMatch(*it);
  }
  matchers().insert(pm);
}

tPDPtr Repository::findParticle(string name) {
  tPDPtr pd;
  string path = name;
  DirectoryAppend(path);
  pd = dynamic_ptr_cast<tPDPtr>(GetPointer(path));
  if ( pd ) return pd;
  for ( ParticleMap::iterator pit = defaultParticles().begin();
	pit != defaultParticles().end(); ++pit )
    if ( pit->second->PDGName() == name ) return pit->second;
  for ( ParticleDataSet::iterator pit = particles().begin();
	pit != particles().end(); ++pit )
    if ( (**pit).PDGName() == name ) return *pit;
  return pd;
}

tPMPtr Repository::findMatcher(string name) {
  for ( MatcherSet::iterator mit = matchers().begin();
	mit != matchers().end(); ++mit )
    if ( name == (**mit).name() ) return *mit;
  return tPMPtr();
}

void Repository::saveRun(string EGname, string name, string filename) {
  EGPtr eg = BaseRepository::GetObject<EGPtr>(EGname);
  EGPtr run = makeRun(eg, name);
  PersistentOStream os(filename, globalLibraries());
  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "Saving event generator '" << name << "'... " << flush;
  os << run;
  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done" << endl;
}

EGPtr Repository::makeRun(tEGPtr eg, string name) {

  // Clone all objects relevant for the EventGenerator. This is
  // the EventGenerator itself, all particles and all particle
  // matchers. 'localObject' is the set of all object refered to by
  // the generator particles and matcher and in the end these are
  // cloned as well.

  // Clone all Particle matchers

  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "Making event generator '" << name << "':" << endl
	   << "Updating all objects... " << flush;

  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done\nCloning matchers and particles... " << flush;

  MatcherSet localMatchers;
  ObjectSet localObjects;
  ObjectSet clonedObjects;
  TranslationMap trans;

  for ( MatcherSet::iterator mit = matchers().begin();
	mit != matchers().end(); ++mit ) {
    PMPtr pm = clone(**mit);
    pm->clear();
    trans[*mit] = pm;
    localMatchers.insert(pm);
    clonedObjects.insert(pm);
    localObjects.insert(*mit);
    addReferences(*mit, localObjects);
  }


  // Clone the particles. But only the ones which should be
  // used. First select the localParticles of the EventGenerator, then
  // add particles from the strategy of the EventGenerator which have
  // not already been selected. Finally add particles from the global
  // default if no default directories has been specified in the
  // strategy which have not already been selected.
  PDVector allParticles;

  for ( ParticleMap::const_iterator pit = eg->localParticles().begin();
 	pit != eg->localParticles().end(); ++pit )
    allParticles.push_back(pit->second);
  if ( eg->strategy() ) {
    tcStrategyPtr strat = eg->strategy();
    for ( ParticleMap::const_iterator pit = strat->particles().begin();
 	  pit != strat->particles().end(); ++pit )
      allParticles.push_back(pit->second);

    vector<string> pdirs;
    if ( eg->strategy()->localParticlesDir().length() )
      pdirs.push_back(eg->strategy()->localParticlesDir());
    pdirs.insert(pdirs.end(), eg->strategy()->defaultParticlesDirs().begin(),
		 eg->strategy()->defaultParticlesDirs().end());

    for ( int i = 0, N = pdirs.size(); i < N; ++i ) {
      string dir = pdirs[i];
      for ( ParticleDataSet::iterator pit = particles().begin();
	    pit != particles().end(); ++pit )
	if ( (**pit).fullName().substr(0, dir.length()) == dir )
	  allParticles.push_back(*pit);
    }
  }

  if ( !eg->strategy() || eg->strategy()->defaultParticlesDirs().empty() )
    for ( ParticleMap::iterator pit = defaultParticles().begin();
	  pit != defaultParticles().end(); ++pit )
      allParticles.push_back(pit->second);

  for ( ParticleDataSet::iterator pit = particles().begin();
	pit != particles().end(); ++pit )
    allParticles.push_back(*pit);

  ParticleMap localParticles;
  set<string> pdgnames;

  for ( PDVector::iterator pit = allParticles.begin();
	pit != allParticles.end(); ++pit ) {
    ParticleMap::iterator it = localParticles.find((**pit).id());
    if ( it == localParticles.end() ) {
      PDPtr pd = clone(**pit);
      trans[*pit] = pd;
      localParticles[pd->id()] = pd;
      clonedObjects.insert(pd);
      localObjects.insert(*pit);
      addReferences(*pit, localObjects);
      if ( pdgnames.find(pd->PDGName()) != pdgnames.end() )
        std::cerr << "Using duplicate PDGName " << pd->PDGName()
                  << " for a new particle.\n This can cause problems and is not "
                  << "recommended.\n If this second particle is a new particle "
                  << "in a BSM Model we recommend you change the name of the particle.\n";
      else
        pdgnames.insert(pd->PDGName());
    } else {
      trans[*pit] = it->second;
    }
  }

  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done\nCloning other objects... " << flush;

  // Clone the OldEventGenerator object to be used:
  localObjects.insert(eg);
  addReferences(eg, localObjects);
  EGPtr egrun = clone(*eg);
  clonedObjects.insert(egrun);
  trans[eg] = egrun;

  for ( ObjectSet::iterator it = localObjects.begin();
	it != localObjects.end(); ++it ) {
    if ( member(trans.map(), *it) ) continue;
    IBPtr ip = clone(**it);
    trans[*it] = ip;
    clonedObjects.insert(ip);
  }

  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done\nRebind references... " << flush;

  IVector defaults;

  trans.translate(inserter(defaults), eg->defaultObjects().begin(),
		  eg->defaultObjects().end());
  if ( eg->strategy() )
    trans.translate(inserter(defaults),
		    eg->strategy()->defaultObjects().begin(),
		    eg->strategy()->defaultObjects().end());

  for ( ObjectSet::iterator it = clonedObjects.begin();
	it != clonedObjects.end(); ++it ) {
    dynamic_cast<Interfaced &>(**it).theGenerator = egrun;
    rebind(**it, trans, defaults);
  }

  // Now, dependencies may have changed, so we do a final round of
  // updates.
  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done\nUpdating cloned objects... " << flush;


  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done\nInitializing... " << flush;

  clonedObjects.erase(egrun);
  egrun->setup(name, clonedObjects, localParticles, localMatchers);

  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done" << endl;

  generators()[name] = egrun;

  return egrun;

}

PDPtr Repository::defaultParticle(PID id) {
  ParticleMap::iterator pit = defaultParticles().find(id);
  return pit == defaultParticles().end()? PDPtr(): pit->second;
}

void Repository::defaultParticle(tPDPtr pdp) {
  if ( pdp ) defaultParticles()[pdp->id()] = pdp;
}

struct ParticleOrdering {
  bool operator()(tcPDPtr p1, tcPDPtr p2) const {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};

struct MatcherOrdering {
  bool operator()(tcPMPtr m1, tcPMPtr m2) const {
    return m1->name() < m2->name() ||
      ( m1->name() == m2->name() && m1->fullName() < m2->fullName() );
  }
};

struct InterfaceOrdering {
  bool operator()(tcIBPtr i1, tcIBPtr i2) const {
    return i1->fullName() < i2->fullName();
  }
};

void Repository::save(string filename) {
  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "saving '" << filename << "'... " << flush;
  PersistentOStream os(filename, globalLibraries());
  set<tcPDPtr,ParticleOrdering>
    part(particles().begin(), particles().end());
  set<tcPMPtr,MatcherOrdering>  match(matchers().begin(), matchers().end());

  os << objects().size();
  for ( ObjectMap::iterator it = objects().begin();
	it != objects().end(); ++it ) os << it->second;
  os << defaultParticles() << part << match << generators()
     << directories() << directoryStack() << globalLibraries() << readDirs();
  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "(" << objects().size() << " objects in " << directories().size()
	   << " directories) done" << endl;
}

string Repository::load(string filename) {
  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "loading '" << filename << "'... " << flush;
  currentFileName() = filename;
  PersistentIStream * is = new PersistentIStream(filename);
  if ( !*is ) {
    delete is;
    // macro is passed in from -D compile flag
    string fullpath = string(THEPEG_PKGLIBDIR) + '/' + filename;
    is = new PersistentIStream(fullpath);
    if ( !*is ) {
      delete is;
      return "Error: Could not find repository '" + filename + "'.";
    }
  }
  *is >> allObjects() >> defaultParticles()
      >> particles() >> matchers() >> generators()
      >> directories() >> directoryStack() >> globalLibraries() >> readDirs();
  delete is;
  objects().clear();
  for ( ObjectSet::iterator it = allObjects().begin();
	it != allObjects().end(); ++it )
    objects()[(**it).fullName()] = *it;

  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "(" << objects().size() << " objects in " << directories().size()
	   << " directories) done\nUpdating... " << flush;
  BaseRepository::resetAll(allObjects());
  BaseRepository::update();

  if ( ThePEG_DEBUG_ITEM(3) )
    clog() << "done" << endl;
  return "";
}

void Repository::stats(ostream & os) {
  os << "number of objects:        " << setw(6) << objects().size() << endl;
  os << "number of objects (all):  " << setw(6) << allObjects().size() << endl;
  os << "number of particles:        " << setw(6) << particles().size() << endl;
  os << "number of matchers:         " << setw(6) << matchers().size() << endl;
}

string Repository::read(string filename, ostream & os) {
  ifstream is;
  string file = filename;
  if ( file[0] == '/' ) {
    if ( ThePEG_DEBUG_LEVEL > 1 ) os << "(= trying to open " << file << " =)" << endl;
    is.open(file.c_str());
  }
  else {
    vector<string> dirs(readDirs().rbegin(), readDirs().rend());
    dirs.push_back(currentReadDirStack().top()); 
    if ( ThePEG_DEBUG_LEVEL > 1 ) {
      os << "(= search path order =)\n(== ";
      std::copy(dirs.rbegin(), dirs.rend(), std::ostream_iterator<string>(os, " ==)\n(== "));
      os << ")" << endl;
    }
    while ( dirs.size() ) {
      string dir = dirs.back();
      if ( dir != "" && dir[dir.length() -1] != '/' ) dir += '/';
      file = dir + filename;
      is.clear();
      if ( ThePEG_DEBUG_LEVEL > 1 ) os << "(= trying to open " << file << " =)" << endl;
      is.open(file.c_str());
      if ( is ) break;
      if ( ThePEG_DEBUG_LEVEL > 1 ) os << "(= no, try next search path =)" << endl;
      dirs.pop_back();
    }
  }
  if ( !is ) {
    return "Error: Could not find input file '" + filename + "'";
  }
  if ( ThePEG_DEBUG_LEVEL > 1 ) os << "(= yes =)" << endl;
  const string dir = StringUtils::dirname(file);
  if ( ThePEG_DEBUG_LEVEL > 1 ) os << "(= pushing <" << dir << "> to stack =)" << endl;
  currentReadDirStack().push(dir);
  try {
    Repository::read(is, os);
    if ( ThePEG_DEBUG_LEVEL > 1 ) os << "(= popping <" << currentReadDirStack().top() << "> from stack =)" << endl;
    currentReadDirStack().pop();
  }
  catch ( ... ) {
    if ( ThePEG_DEBUG_LEVEL > 1 ) os << "(= popping <" << currentReadDirStack().top() << "> from stack =)" << endl;
    currentReadDirStack().pop();
    throw;
  }
  return "";
}

string Repository::
modifyEventGenerator(EventGenerator & eg, string filename, 
		     ostream & os, bool initOnly) {
  ObjectSet objs = eg.objects();
  objs.insert(&eg);
  for ( ObjectSet::iterator it = objs.begin(); it != objs.end(); ++it ) {
    string name = (**it).fullName();
    if ( name.rfind('/') != string::npos )
      CreateDirectory(name.substr(0, name.rfind('/') + 1));
    objects()[name] = *it;
    allObjects().insert(*it);
  }
  
  string msg = read(filename, os);

  if ( !msg.empty() )
    return msg;
 
  for_each(objs, mem_fn(&InterfacedBase::reset));
  eg.initialize(initOnly);

  if ( !generators().empty() )
    msg += "Warning: new generators were initialized while modifying "
      + eg.fullName() + ".\n";

  return msg;
}

void Repository::resetEventGenerator(EventGenerator & eg) {

  ObjectSet objs = eg.objects();
  objs.insert(&eg);
  for ( ObjectSet::iterator it = objs.begin(); it != objs.end(); ++it ) {
    string name = (**it).fullName();
    if ( name.rfind('/') != string::npos )
      CreateDirectory(name.substr(0, name.rfind('/') + 1));
    objects()[name] = *it;
    allObjects().insert(*it);
  }
  
  for_each(objs, mem_fn(&InterfacedBase::reset));
  eg.initialize(true);

}

void Repository::execAndCheckReply(string line, ostream & os) {
  string reply = exec(line, os);
  if ( reply.size() ) 
    os << reply;
  if ( reply.size() && reply[reply.size()-1] != '\n' ) 
    os << endl;
  if ( exitOnError() && reply.size() >= 7 
       && reply.substr(0, 7) == "Error: " )
    exit(exitOnError());
}

void Repository::read(istream & is, ostream & os, string prompt) {
#ifdef HAVE_LIBREADLINE
  if ( &is == &std::cin ) {
    char * line_read = 0;
    do {
      if ( line_read ) {
	free(line_read);
	line_read = 0;
      }
      
      line_read = readline(prompt.c_str());
      
      if ( line_read && *line_read ) {
	string line = line_read;
	while ( !line.empty() && line[line.size() - 1] == '\\' ) {
	  line[line.size() - 1] = ' ';
	  char * cont_read = readline("... ");
	  if ( cont_read ) {
	    line += cont_read;
	    free(cont_read);
	  }
	}
	if ( prompt.empty() && ThePEG_DEBUG_LEVEL > 0 )
	  os << "(" << line << ")" << endl;
#ifdef HAVE_READLINE_HISTORY
	add_history(line.c_str());
#endif // HAVE_READLINE_HISTORY
	execAndCheckReply(line, os);
      }
    }
    while ( line_read );
  }
  else {
#endif // HAVE_LIBREADLINE
    string line;
    if ( prompt.size() ) os << prompt;
    while ( getline(is, line) ) {
      while ( !line.empty() && line[line.size() - 1] == '\\' ) {
	line[line.size() - 1] = ' ';
	string cont;
	if ( prompt.size() ) os << "... ";
	getline(is, cont);
	line += cont;
      }
      if ( prompt.empty() && ThePEG_DEBUG_LEVEL > 0 )
	os << "(" << line << ")" << endl;
      execAndCheckReply(line, os);
      if ( prompt.size() ) os << prompt;
    }
#ifdef HAVE_LIBREADLINE
  }
#endif
  if ( prompt.size() ) os << endl;
}


string Repository::copyParticle(tPDPtr p, string newname) {
  DirectoryAppend(newname);
  
  string newdir = newname.substr(0, newname.rfind('/')+1);
  newname =newname.substr(newname.rfind('/')+1);
  if ( newname.empty() ) newname = p->name();
  if ( GetPointer(newdir + newname) )
    return "Error: Cannot create particle " + newdir + newname +
      ". Object already exists.";
  if ( p->CC() && GetPointer(newdir + p->CC()->name()) )
    return "Error: Cannot create anti-particle " + newdir + newname +
      ". Object already exists.";
  PDPtr pd = p->pdclone();
  Register(pd, newdir + newname);
  pd->theDecaySelector.clear();
  pd->theDecayModes.clear();
  pd->isStable = true;
  if ( p->CC() ) {
    PDPtr apd = p->CC()->pdclone();
    Register(apd, newdir + apd->name());
    apd->theDecaySelector.clear();
    apd->theDecayModes.clear();
    apd->isStable = true;
    pd->theAntiPartner = apd;
    apd->theAntiPartner = pd;
    pd->syncAnti = p->syncAnti;
    apd->syncAnti = p->CC()->syncAnti;
  }
  HoldFlag<> dosync(pd->syncAnti, true);
  for ( DecaySet::const_iterator it = p->theDecayModes.begin();
	it != p->theDecayModes.end(); ++it )
    pd->addDecayMode(*it);
  return "";
}

void Repository::remove(tIBPtr ip) {
  ObjectMap::iterator it = objects().find(ip->fullName());
  if ( it == objects().end() || ip != it->second ) return;
  objects().erase(it);
  allObjects().erase(ip);
  if ( dynamic_ptr_cast<tPDPtr>(ip) ) {
    particles().erase(dynamic_ptr_cast<tPDPtr>(ip));
    defaultParticles().erase(dynamic_ptr_cast<tPDPtr>(ip)->id());
  }
  if ( dynamic_ptr_cast<tPMPtr>(ip) )
    matchers().erase(dynamic_ptr_cast<tPMPtr>(ip));
}

string Repository::remove(const ObjectSet & rmset) {
  ObjectSet refset;
  for ( ObjectMap::const_iterator i = objects().begin();
	i != objects().end(); ++i ) {
    if ( member(rmset, i->second) ) continue;
    IVector ov = DirectReferences(i->second);
    for ( int j = 0, M = ov.size(); j < M; ++j )
      if ( member(rmset, ov[j]) ) {
	refset.insert(i->second);
	break;
      }
  }
  if ( refset.empty() ) {
    for ( ObjectSet::iterator oi = rmset.begin(); oi != rmset.end(); ++oi )
      remove(*oi);
    return "";
  }
  string ret = "Error: cannot remove the objects because the following "
    "objects refers to some of them:\n";
  for ( ObjectSet::iterator oi = refset.begin(); oi != refset.end(); ++oi )
    ret += (**oi).fullName() + "\n";
  return ret;
}
   
string Repository::exec(string command, ostream & os) {
  string cpcmd = command;
  try {
    string verb = StringUtils::car(command);
    command = StringUtils::cdr(command);
    if ( verb == "help" ) {
      help(command, os);
      return "";
    }
    if ( verb == "rm" ) {
      ObjectSet rmset;
      while ( !command.empty() ) {
	string name = StringUtils::car(command);
	DirectoryAppend(name);
	IBPtr obj = GetPointer(name);
	if ( !obj ) return "Error: Could not find object named " + name;
	rmset.insert(obj);
	command = StringUtils::cdr(command);
      }
      return remove(rmset);
    }
    if ( verb == "rmdir" || verb == "rrmdir" ) {
      string dir = StringUtils::car(command);
      DirectoryAppend(dir);
      if ( dir[dir.size() - 1] != '/' ) dir += '/';
      if ( !member(directories(), dir) )
	return verb == "rmdir"? "Error: No such directory.": "";
      IVector ov = SearchDirectory(dir);
      if ( ov.size() && verb == "rmdir" )
	return "Error: Cannot remove a non-empty directory. "
	  "(Use rrmdir do remove all object and subdirectories.)";
      ObjectSet rmset(ov.begin(), ov.end());
      string ret = remove(rmset);
      if ( !ret.empty() ) return ret;
      StringVector dirs(directories().begin(), directories().end());
      for ( int i = 0, N = dirs.size(); i < N; ++ i )
	if ( dirs[i].substr(0, dir.size()) == dir )
	  directories().erase(dirs[i]);
      for ( int i = 0, N = directoryStack().size(); i < N; ++i )
	if ( directoryStack()[i].substr(0, dir.size()) == dir )
	  directoryStack()[i] = '/';
      return "";
    }
    if ( verb == "cp" ) {
      string name = StringUtils::car(command);
      DirectoryAppend(name);
      tPDPtr p = dynamic_ptr_cast<tPDPtr>(GetPointer(name));
      if ( p ) return copyParticle(p, StringUtils::cdr(command));
      return BaseRepository::exec(cpcmd, os);
    }
    if ( verb == "setup" ) {
      string name = StringUtils::car(command);
      DirectoryAppend(name);
      IBPtr obj = GetPointer(name);
      if ( !obj ) return "Error: Could not find object named " + name;
      istringstream is(StringUtils::cdr(command));
      readSetup(obj, is);
      // A particle may have been registered before but under the wrong id().
      PDPtr pd = dynamic_ptr_cast<PDPtr>(obj);
      if(pd) registerParticle(pd);
      return "";
    }
    if ( verb == "decaymode" ) {
      string tag = StringUtils::car(command);
      DMPtr dm = DecayMode::constructDecayMode(tag);
      if ( !dm ) return "Error: Could not create decay mode from the tag " +
		   StringUtils::car(command);
      istringstream is(StringUtils::cdr(command));
      readSetup(dm, is);
      if ( !dm->CC() ) return "";

      if ( dm->CC()->parent()->synchronized() ) {
	dm->CC()->synchronize();
	return "";
      }

      if ( !dm->CC()->decayer() )
	return FindInterface(dm, "Decayer")->
	  exec(*dm->CC(), "set", dm->decayer()->fullName());
      return "";
    }
    if ( verb == "makeanti" ) {
      string name = StringUtils::car(command);
      DirectoryAppend(name);
      tPDPtr p = dynamic_ptr_cast<tPDPtr>(GetPointer(name));
      if ( !p ) return "Error: No particle named " + name;
      name = StringUtils::car(StringUtils::cdr(command));
      DirectoryAppend(name);
      tPDPtr ap = dynamic_ptr_cast<tPDPtr>(GetPointer(name));
      if ( !ap ) return "Error: No particle named " + name;
      ParticleData::antiSetup(PDPair(p, ap));
      return "";
    }
    if ( verb == "read" ) {
      // remember directory we're in
      string cwd = directoryStack().back();
      string filename = StringUtils::car(command);
      string msg = read(filename, os);
      // Return to the original directory, so that
      // calling 'read' in an input file will not change the 
      // repository directory you're in
      ChangeDirectory(cwd);
      return msg;
    }
    if ( verb == "load" ) {
      return load(StringUtils::car(command));
    }      
    if ( verb == "save" ) {
      save(StringUtils::car(command));
      return "";
    }
    if ( verb == "lsruns" ) {
      string ret;
      for ( GeneratorMap::iterator ieg = generators().begin();
	    ieg != generators().end(); ++ieg ) ret += ieg->first + "\n";
      return ret;
    }
    if ( verb == "makerun" ) {
      string runname = StringUtils::car(command);
      string generator = StringUtils::car(StringUtils::cdr(command));
      DirectoryAppend(generator);
      EGPtr eg = BaseRepository::GetObject<EGPtr>(generator);
      makeRun(eg, runname);
      return "";
    }
    if ( verb == "rmrun" ) {
      string runname = StringUtils::car(command);
      generators().erase(runname);
      return "";
    }
    if ( verb == "saverun" || verb == "saverunfile" || verb == "run" ) {
      string runname = StringUtils::car(command);
      string generator = StringUtils::car(StringUtils::cdr(command));
      DirectoryAppend(generator);
      GeneratorMap::iterator ieg = generators().find(runname);
      EGPtr eg;
      if ( ieg == generators().end() ) {
	eg = BaseRepository::GetObject<EGPtr>(generator);
	eg = makeRun(eg, runname);
      } else
	eg = ieg->second;
      if ( !eg )
	return "Error: Could not create/find run named'" + runname + "'.";
      if ( verb == "run" ) 
	eg->go();
      else if ( verb == "saverunfile" ) {
	string file = generator;
	PersistentOStream os(file, globalLibraries());
	os << eg;
	if ( !os ) return "Save failed! (I/O error)";
      } else {
	string file = eg->filename() + ".run";
	PersistentOStream os(file, globalLibraries());
	os << eg;
	if ( !os ) return "Save failed! (I/O error)";
      }
      return "";
    }
    if ( verb == "removerun" ) {
      string runname = StringUtils::car(command);
      GeneratorMap::iterator ieg = generators().find(runname);
      if ( ieg != generators().end() ) {
	generators().erase(ieg);
	return "";
      } else
	return "Error: No run named '" + runname + "' available.";
    }
    if ( verb == "create" ) {
      string className = StringUtils::car(command);
      command = StringUtils::cdr(command);
      string name = StringUtils::car(command);
      const ClassDescriptionBase * db = DescriptionList::find(className);
      command = StringUtils::cdr(command);
      while ( !db && command.length() ) {
	string library = StringUtils::car(command);
	command = StringUtils::cdr(command);
	DynamicLoader::load(library);
	db = DescriptionList::find(className);
      }
      if ( !db ) {
	string msg = "Error: " + className + ": No such class found.";
	if ( !DynamicLoader::lastErrorMessage.empty() )
	  msg += "\nerror message from dynamic loader:\n" +
	    DynamicLoader::lastErrorMessage;
	return msg;
      }
      IBPtr obj = dynamic_ptr_cast<IBPtr>(db->create());
      if ( !obj ) return "Error: Could not create object of this class class.";
      if ( name.empty() ) return "Error: No name specified.";
      Register(obj, name);
      return "";
    }
    if ( verb == "defaultparticle" ) {
      while ( !command.empty() ) {
	string name = StringUtils::car(command);
	DirectoryAppend(name);
	tPDPtr p = dynamic_ptr_cast<tPDPtr>(GetPointer(name));
	if ( !p ) return "Error: No particle named " + name;
	defaultParticle(p);
	command = StringUtils::cdr(command);
      }
      return "";
    }
    if ( verb == "EXITONERROR" ) {
      exitOnError() = 1;
      return "";
    }
  }
  catch (const Exception & e) {
    e.handle();
    return "Error: " + e.message();
  }

  return BaseRepository::exec(cpcmd, os);
}

void Repository::help(string cmd, ostream & os) {
 
  cmd = StringUtils::car(cmd);

  if ( cmd == "cd" )
    os << "Usage: cd <directory>" << endl
       << "Set the current directory to <directory>." << endl;
  else if ( cmd == "mkdir" )
    os << "Usage: mkdir <path-name>" << endl
       << "Create a new directory called with the given path name." << endl;
  else if ( cmd == "rmdir" )
    os << "Usage: rmdir <directory>" << endl
       << "Remove an empty directory." << endl;
  else if ( cmd == "rrmdir" )
    os << "Usage: rrmdir <directory>" << endl
       << "Remove a directory and everything that is in it recursively." << endl
       << "Will only succeed if no other objects refers to the ones to "
       << "be deleted." << endl;
  else if ( cmd == "cp" )
    os << "Usage: cp <object> <path-name>" << endl
       << "Copy the given object to a new object with the given name." << endl;
  else if ( cmd == "setup" )
    os << "Usage: setup <object> <arguments> ..." << endl
       << "Tell a given object to read information given by the arguments."
       << endl;
  else if ( cmd == "decaymode" )
    os << "Usage: decaymode <tag> <branching fraction> <on|off> <decayer-object>"
       << endl
       << "Construct a decay mode from the given decay tag. The resulting "
       << "object will be inserted in the directory with the same path as "
       << "the decaying particle object. The given brancing fraction will "
       << "be set as well as the given decayer object. If the mode should "
       << "be switched on by default 1(on) should be specified (otherwise "
       << "0(off))." << endl;
  else if ( cmd == "makeanti" )
    os << "Usage: makeanti <particle-object> <particle-object>" << endl
       << "Indicate that the two given particle objects are eachothers "
       << "anti-partnets." << endl;
  else if ( cmd == "read" )
    os << "Usage: read <file-name>" << endl
       << "Read more commands from the given file. The file name can be "
       << "given relative to the current directory in the shell, or "
       << "relative to standard directories, or as an absolute path." << endl;
  else if ( cmd == "load" )
    os << "Usage: load <repository-file-name>" << endl
       << "Discard everything in the reopsitory and read in a completely "
       << "new repository from the given file." << endl;
  else if ( cmd == "save" )
    os << "Usage: save <file-name>" << endl
       << "Save the complete repository to the given file." << endl;
  else if ( cmd == "lsruns" )
    os << "Usage: lsruns" << endl
       << "List the run names of all initialized event generators." << endl;
  else if ( cmd == "makerun" )
    os << "Usage: makerun <run-name> <event-generator-object>" << endl
       << "Initialize the given event generator and assign a run name." << endl;
  else if ( cmd == "rmrun" )
    os << "Usage: rmrun <run-name>" << endl
       << "Remove the initialized event generator given by the run name."
       << endl;
  else if ( cmd == "saverun" )
    os << "Usage: saverun <run-name> <event-generator-object>" << endl
       << "Initialize the given event generator and assign a run name "
       << "and save it to a file named <run-name>.run" << endl;
  else if ( cmd == "run" )
    os << "Usage: run <run-name>" << endl
       << "Run the initialized event generator given b the run name." << endl;
  else if ( cmd == "create" )
    os << "Usage: create <class-name> <name> {<dynamic-library>}" << endl
       << "Create an object of the given class and assign the given name. "
       << "Optionally supply a dynamically loaded library where the class "
       << "is included." << endl;
  else if ( cmd == "pushd" )
    os << "Usage: pushd <directory>" << endl
       << "Set the current directory to <directory>, but keep the previous "
       << "working directory on the directory stack." << endl;
  else if ( cmd == "popd" )
    os << "Usage: popd" << endl
       << "Leave the current working directory and set the current "
       << "directory to the previous one on the directory stack." << endl;
  else if ( cmd == "pwd" )
    os << "Usage: pwd" << endl
       << "Print the current working directory." << endl;
  else if ( cmd == "dirs" )
    os << "Usage: dirs" << endl
       << " Print the contents of the directory stack." << endl;
  else if ( cmd == "mv" )
    os << "Usage: mv  <object> <path-name>" << endl
       << "Rename the given object to a new path name." << endl;
  else if ( cmd == "ls" )
    os << "Usage: ls {<directory>}" << endl
       << "List the objects and subdirectories in the current or given "
       << "directory." << endl;
  else if ( cmd == "library" )
    os << "Usage: library <dynamic-library>" << endl
       << "Make new classes available to the repository by dynamically "
       << "linking the given library." << endl;
  else if ( cmd == "globallibrary" )
    os << "Usage: globallibrary <dynamic-library>" << endl
       << "Make new classes available to the repository by dynamically "
       << "linking the given library. If this repository is saved and read "
       << "in again, this library will be linked in from the beginning." << endl;
  else if ( cmd == "rmgloballibrary" )
    os << "Usage: rmgloballibrary <dynamic-library>" << endl
       << "Remove a dynamic library previously added with globallibrary."
       << endl;
  else if ( cmd == "appendpath" )
    os << "Usage: appendpath <unix-directory>" << endl
       << "Add a search path for dynamic libraries to the end of the "
       << "search list." << endl;
  else if ( cmd == "lspaths" )
    os << "Usage: lspaths" << endl
       << "List search paths for dynamic libraries." << endl;
  else if ( cmd == "prependpath" )
    os << "Usage: prependpath <unix-directory>" << endl
       << "Add a search path for dynamic libraries to the beginning of the "
       << "search list." << endl;
  else if ( cmd == "doxygendump" )
    os << "Usage: doxygendump <namespace> <filename>" << endl
       << "Extract doxygen documentation of all loaded classes in the "
       << "given name space and weite it to a file.." << endl;
  else if ( cmd == "mset" || cmd == "minsert" || cmd == "mdo" )
    os << "Usage: " << cmd << " <directory> <class> <interface> <value>" << endl
       << "Recursively find in the given directory all objects of the "
       << "given class and call '" << cmd.substr(1)
       << "' with the given value for the given interface." << endl;
  else if ( cmd == "msetdef" || cmd == "mget" || cmd == "mdef" ||
	    cmd == "mmin" || cmd == "mmax" || cmd == "merase" )
    os << "Usage: " << cmd << " <directory> <class> <interface>" << endl
       << "Recursively find in the given directory all objects of the given "
       << "class and call '" << cmd.substr(1)
       << "' for the given interface." << endl;
  else if ( cmd == "set" )
    os << "Usage: set <object>:<interface> <value>" << endl
       << "Set the interface for the given object to the given value." << endl;
  else if ( cmd == "setdef" )
    os << "Usage: setdef <object>:<interface>" << endl
       << "Set the interface for the given object to its default value." << endl;
  else if ( cmd == "insert" )
    os << "Usage: insert <object>:<interface> <value>" << endl
       << "Insert a value in the vector interface of the given object." << endl;
  else if ( cmd == "erase" )
    os << "Usage: erase <object>:<interface>" << endl
       << "Erase a value from the vector interface of the given object." << endl;
  else if ( cmd == "do" )
    os << "Usage: do <object>:<command-interface> <arguments>" << endl
       << "Call the command interface of the given object with the "
       << "given arguments." << endl;
  else if ( cmd == "get" )
    os << "Usage: get <object>:<interface>" << endl
       << "Print the value of the interface of the given object." << endl;
  else if ( cmd == "def" )
    os << "Usage: def <object>:<interface>" << endl
       << "Print the default value of the interface of the given object."
       << endl;
  else if ( cmd == "min" )
    os << "Usage: min <object>:<interface>" << endl
       << "Print the minimum value of the interface of the given object."
       << endl;
  else if ( cmd == "max" )
    os << "Usage: max <object>:<interface>" << endl
       << "Print the maximum value of the interface of the given object."
       << endl;
  else if ( cmd == "describe" )
    os << "Usage: describe <object>{:<interface>}" << endl
       << "Describe the given object or an interface of the object." << endl;
  else if ( cmd == "lsclass" )
    os << "Usage: lsclass" << endl
       << "List all classes available in the repository." << endl;
  else if ( cmd == "all" ) {
    os << "Available commands:"
       << endl
       << "* cd, mkdir, rmdir, rrmdir, pwd, cp, mv, rm, pushd, popd, dirs, ls:\n"
       << "  Manipulate the repository structure. Analogous to unix "
       << "shell commands."
       << endl
       << "* create, setup, decaymode makeanti:\n"
       << "  Create or setup an object."
       << endl
       << "* set, get, insert, erase, do, detdef, def, min, max, describe\n"
       << "  mset, minsert, mdo, msetdef, mdef, mmin, mmax, merase:\n"
       << "  Manipulate interfaces to objects."
       << endl
       << "* makerun, saverun, run, lsruns, rmrun:\n"
       << "  Create and handle initialized event genrators which can be run."
       << endl
       << "* read, load, library globallibrary, rmgloballibrary,\n"
       << "  appendpath, prependpath, lspaths, doxygendump:\n"
       << "  Handle files external files and libraries."
       << endl;
    os << "Do 'help syntax' for help on syntax." << endl
       << "Do 'help <command>' for help on a particular command." << endl;
  }
  else if ( cmd == "syntax" )
    os << "* <directory> = '/' | <name> | <directory>/<name>" << endl
       << "  <object> = <name> | <directory>/<name> | <object>:<ref-interface>\n"
       << "  Analogous to a unix file structure, an object can be "
       << "specified with an\n  absolute path or a path relative to "
       << "the current directory." << endl
       << "* <interface> = <interface-name>|<interface-name>[<index>]" << endl
       << "  An interface can be a parameter (floating point, integer or "
       << "string),\n  a switch (integer, possibly named), a reference to "
       << "another object in the\n  repository or a command which takes "
       << "an arbitrary string as argument.\n  There are also vector interfaces "
       << "of parameters and references for which\n  an index must be supplied."
       << endl;
  else {
    if ( !cmd.empty() ) os << "No command '" << cmd << "' found." << endl;
    os << "Common commands:" << endl
       << "* cd, mkdir, rmdir, pwd, cp, mv, rm:\n"
       << "  Manipulate the repository structure. Analogous to unix "
       << "shell commands." << endl
       << "* create, setup:\n"
       << " Create an object." << endl
       << "set, get, insert, erase, do:\n"
       << " Manipulate interfaces to objects." << endl
       << "* makerun, saverun, run, lsruns:\n"
       << " Create and handle initialized event genrators which can be run."
       << endl;
    os << "Do 'help all' for a complete list of commands." << endl
       << "Do 'help syntax' for help on syntax." << endl
       << "Do 'help <command>' for help on a particular command." << endl;
  }

}

Repository::Repository() {
  ++ninstances;
}

Repository::~Repository() {
  --ninstances;
  if ( ninstances <= 0 ) {
    generators().clear();
  }
}

int Repository::ninstances = 0;

namespace {
static string version_ =
#include "versionstamp.inc"
"";
}

string Repository::version() {
  return ::version_;
}

string Repository::banner() {
  const auto now    = std::chrono::system_clock::now();
  const auto now_c  = std::chrono::system_clock::to_time_t(now);
  string time = ">>>> " ;
  time += StringUtils::stripws(string(std::ctime(&now_c))) + ' ';
  time += string(max(0,74 - int(time.size())), ' ');
  time += "<<<<";

  string line = ">>>> Toolkit for HEP Event Generation - "
    + Repository::version() + ' ';
  line += string(max(0,78 - int(line.size())), '<');

  string block = string(78, '>') + '\n' 
                 + line + '\n' 
                 + time + '\n'
                 + string(78, '<') + '\n';
  return block;
}

