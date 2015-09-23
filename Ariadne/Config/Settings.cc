// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Settings class.
//

#include "Settings.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/Repository.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

Settings::~Settings() {}

IBPtr Settings::clone() const {
  return new_ptr(*this);
}

IBPtr Settings::fullclone() const {
  return new_ptr(*this);
}

bool Settings::preInitialize() const {
  return true;
}

void Settings::doupdate() {
  HandlerBase::doupdate();
  // First update base class.
  bool redo = touched();
  // redo if touched.
//  UpdateChecker::check(aDependentMember, redo);
  // Update referenced objects on which this depends redo is set to true
  // if the dependent object is touched.
//  for_each(ContainerOfDependencies, UpdateChecker(redo));
  // Update a container of references.
//  for_each(MapOfDependencies, UpdateMapChecker(redo));
  // Update a map of references.
  if ( !redo ) return;
  // return if nothing has been touched. Otherwise do the actual update.
//  touch()
  // Touch if anything has changed.
}

void Settings::doinit() {
  HandlerBase::doinit();
  for ( int i = 0, N = setObjects.size(); i < N; ++i )
    generator()->preinitInterface(dynamic_ptr_cast<IPtr>(setObjects[i]),
				  setInterfaces[i], "set", setValues[i]);
}

void Settings::dofinish() {
  HandlerBase::dofinish();
}

void Settings::doinitrun() {
  HandlerBase::doinitrun();
  if ( !verb ) return;
  generator()->log() << "> Settings object " << fullName()
		     << " has set the following interfaces:" << endl;
  for ( int i = 0, N = setObjects.size(); i < N; ++i )
    generator()->log() << "set " << setObjects[i]->fullName() << ":" << setInterfaces[i]
		     << " " << setValues[i] << endl;
  generator()->log() << endl;
}

template <typename Cont>
void translate(const TranslationMap & trans, Cont & c) {
  Cont old;
  old.swap(c);
  trans.translate(inserter(c), old.begin(), old.end());
}

void Settings::rebind(const TranslationMap & trans) {
  // dummy = trans.translate(dummy);
  HandlerBase::rebind(trans);
  translate(trans, setObjects);
  translate(trans, setRefs);
}

IVector Settings::getReferences() {
  IVector ret = HandlerBase::getReferences();
  ret.insert(ret.end(), setObjects.begin(), setObjects.end());
  ret.insert(ret.end(), setRefs.begin(), setRefs.end());
  return ret;
}

string Settings::setFunction(string cmd) {
  string noun = StringUtils::car(cmd);
  IBPtr ip = Repository::getObjectFromNoun(noun);
  const InterfaceBase * ifb =
    Repository::FindInterface(ip, Repository::getInterfaceFromNoun(noun));
  if ( !ip ) return "Object not found!";
  if ( !ifb ) return "Interface not found!";
  string arg = StringUtils::cdr(cmd);
  if ( dynamic_cast<const RefInterfaceBase *>(ifb) ) {
    IBPtr ipr = Repository::getObjectFromNoun(arg);
    if ( ipr ) {
      arg = ipr->fullName();
      setRefs.push_back(ipr);
    }
  }
  string arguments = Repository::getPosArgFromNoun(noun) + " "
    + arg;
  setObjects.push_back(ip);
  setInterfaces.push_back(ifb->name());
  setValues.push_back(arguments);
  return"";
}

void Settings::persistentOutput(PersistentOStream & os) const {
  os << setObjects << setInterfaces << setValues << setRefs << verb;
}

void Settings::persistentInput(PersistentIStream & is, int) {
  is >> setObjects >> setInterfaces >> setValues >> setRefs >> verb;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Settings,HandlerBase>
  describeThePEGSettings("ThePEG::Settings", "libAriadne5.so");

void Settings::Init() {

  static ClassDocumentation<Settings> documentation
    ("There is no documentation for the Settings class");


  static Command<Settings> interfaceset
    ("set",
     "Set an interface for an object at initialization.",
     &Settings::setFunction, true);

  static Switch<Settings,bool> interfaceVerbosity
    ("Verbosity",
     "Optionally write in the log-file what settings has been done.",
     &Settings::verb, true, true, false);
  static SwitchOption interfaceVerbosityVerbose
    (interfaceVerbosity,
     "Verbose",
     "Write setting to log-file.",
     true);
  static SwitchOption interfaceVerbosityQuiet
    (interfaceVerbosity,
     "Quiet",
     "Do not write in the log-file.",
     false);

}

