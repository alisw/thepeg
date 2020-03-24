// -*- C++ -*-
//
// BaseRepository.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaseRepository class.
//

// macro is passed in from -D compile flag
#ifndef THEPEG_PKGDATADIR
#error Makefile.am needs to define THEPEG_PKGDATADIR
#endif

#include "BaseRepository.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/InterfaceBase.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/Utilities/DescriptionList.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Utilities/TypeInfo.h"
#include "ThePEG/Utilities/DynamicLoader.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "BaseRepository.tcc"
#endif

using namespace ThePEG;

ostream *& BaseRepository::coutp() {
  static ostream * theCout = &std::cout;
  return theCout;
}

ostream *& BaseRepository::cerrp() {
  static ostream * theCerr = &std::cerr;
  return theCerr;
}

ostream *& BaseRepository::clogp() {
  static ostream * theClog = &std::clog;
  return theClog;
}

bool & BaseRepository::updating() {
  static bool theBool = false;
  return theBool;
}

ObjectMap & BaseRepository::objects() {
  static ObjectMap theObjectMap;
  return theObjectMap;
}

ObjectSet & BaseRepository::allObjects() {
  static ObjectSet theObjectSet;
  return theObjectSet;
}

BaseRepository::TypeInterfaceMap & BaseRepository::interfaces() {
  static TypeInterfaceMap theInterfaceMap;
  return theInterfaceMap;
}

BaseRepository::TypeDocumentationMap & BaseRepository::documentations() {
  static TypeDocumentationMap theDocumentationMap;
  return theDocumentationMap;
}

BaseRepository::DirectorySet & BaseRepository::directories() {
  static DirectorySet theDirectories = {"/"};
  return theDirectories;
}

vector<string> & BaseRepository::globalLibraries() {
  static vector<string> theGlobalLibraries;
  return theGlobalLibraries;
}

stack<string> & BaseRepository::currentReadDirStack() {
  static stack<string> theCurrentReadDirStack;
  if ( theCurrentReadDirStack.empty() ) theCurrentReadDirStack.push("");
  return theCurrentReadDirStack;
}

vector<string> & BaseRepository::readDirs() {
  // macro is passed in from -D compile flag
  static vector<string>
    theReadDirs(1, THEPEG_PKGDATADIR);
  return theReadDirs;
}

const vector<string> & BaseRepository::listReadDirs() {
  return BaseRepository::readDirs();
}

void BaseRepository::prependReadDir(string dir) {
  readDirs().insert(readDirs().begin(), dir);
}

void BaseRepository::prependReadDir(const std::vector<std::string>& dirs) {
  readDirs().insert(readDirs().begin(), dirs.begin(), dirs.end());
}

void BaseRepository::appendReadDir(string dir) {
  readDirs().push_back(dir);
}

void BaseRepository::appendReadDir(const std::vector<std::string>& dirs) {
  readDirs().insert(readDirs().end(), dirs.begin(), dirs.end());
}

BaseRepository::StringVector & BaseRepository::directoryStack() {
  static StringVector theDirectoryStack(1, "/");
  return theDirectoryStack;
}

void BaseRepository::Register(const InterfaceBase & ib, const type_info & i) {
  const ClassDescriptionBase * db = DescriptionList::find(i);
  if ( db ) interfaces()[db].insert(&ib);
}

void BaseRepository::
Register(const ClassDocumentationBase & cd, const type_info & i) {
  const ClassDescriptionBase * db = DescriptionList::find(i);
  if ( db ) documentations()[db] = &cd;
}

void BaseRepository::Register(IBPtr ip, string newName) {
  DirectoryAppend(newName);
  ip->name(newName);
  Register(ip);
}

void BaseRepository::Register(IBPtr ip) {
  if ( !ip || member(allObjects(), ip) ) return;
  if ( member(objects(), ip->fullName()) )
    throw RepoNameExistsException(ip->fullName());
  objects()[ip->fullName()] = ip;
  allObjects().insert(ip);
  ip->clear();
  ip->update();
  ip->touch();
}

void BaseRepository::DirectoryAppend(string & name) {
  if ( name == "." ) name = directoryStack().back();
  if ( name[0] != '/' ) name = directoryStack().back() + name;
}  

void BaseRepository::CreateDirectory(string name) {
  DirectoryAppend(name);
  if ( name[name.size()-1] != '/' ) name += "/";
  if ( member(directories(), name) ) return;
  directories().insert(name);
  name = name.substr(0, name.size() - 1);
  name = name.substr(0, name.rfind('/'));
  if ( name.size() ) CreateDirectory(name);
}

void BaseRepository::CheckObjectDirectory(string name) {
  if ( name[name.size() - 1] != '/' )
    name = name.substr(0, name.rfind('/') + 1);
  CheckDirectory(name);
}

void BaseRepository::CheckDirectory(string name) {
  DirectoryAppend(name);
  if ( name[name.size()-1] != '/' ) name += "/";  
  if ( member(directories(), name) ) return;
  throw RepositoryNoDirectory(name);
}

void BaseRepository::ChangeDirectory(string name) {
  DirectoryAppend(name);
  if ( name[name.size()-1] != '/' ) name += "/";  
  if ( member(directories(), name) ) {
    directoryStack().back() = name;
    return;
  }
  throw RepositoryNoDirectory(name);
}

void BaseRepository::PushDirectory(string name) {
  DirectoryAppend(name);
  if ( name[name.size()-1] != '/' ) name += "/";  
  if ( member(directories(), name) ) {
    directoryStack().push_back(name);
    return;
  }
  throw RepositoryNoDirectory(name);
}

void BaseRepository::PopDirectory() {
  if ( directoryStack().size() > 1 ) directoryStack().pop_back();
}

IBPtr BaseRepository::GetPointer(string name) {
  ObjectMap::iterator it = objects().find(name);
  return it == objects().end()? IBPtr(): it->second;
}

IVector BaseRepository::SearchDirectory(string name, string className) {
  IVector ret;
  DirectoryAppend(name);
  const ClassDescriptionBase * cdb = 0;
  if ( className.size() ) {
    cdb = DescriptionList::find(className);
    if ( !cdb ) return ret;
  }
  if ( name[name.size()-1] != '/' ) name += "/";
  string::size_type size = name.size();
  for ( ObjectMap::const_iterator i = objects().begin();
	i != objects().end(); ++i ) {
    const auto & tmp=(*(i->second));
    if ( cdb && !DescriptionList::find(typeid(tmp))->isA(*cdb) )
      continue;
    if ( i->first.substr(0, size) == name ) ret.push_back(i->second);
  }
  return ret;
}

IVector BaseRepository::GetObjectsReferringTo(IBPtr obj) {
  IVector ret;
  for ( ObjectMap::const_iterator i = objects().begin();
	i != objects().end(); ++i ) {
    if ( obj == i->second ) continue;
    IVector ov = DirectReferences(i->second);
    if ( member(ov, obj) ) ret.push_back(i->second);
  }
  return ret;
}

IVector BaseRepository::DirectReferences(IBPtr obj) {
  IVector ov = obj->getReferences();
  const auto & tmp=*obj;
  InterfaceMap interfaceMap = getInterfaces(typeid(tmp));
  for ( InterfaceMap::iterator iit = interfaceMap.begin();
	iit != interfaceMap.end(); ++iit ) {
    IVector ovi = iit->second->getReferences(*obj);
    ov.insert(ov.end(), ovi.begin(), ovi.end());
  }
  return ov;
}

void BaseRepository::
addReferences(tIBPtr obj, ObjectSet & refs) {
  if ( !obj ) return;
  refs.insert(obj);
  IVector ov = obj->getReferences();
  for ( IVector::const_iterator it = ov.begin(); it != ov.end(); ++it )
    if ( !member(refs, *it) ) addReferences(*it, refs);
  const auto & tmp=*obj;
  InterfaceMap interfaceMap = getInterfaces(typeid(tmp));
  for ( InterfaceMap::iterator iit = interfaceMap.begin();
	iit != interfaceMap.end(); ++iit ) {
    IVector ov = iit->second->getReferences(*obj);
    for ( IVector::const_iterator it = ov.begin(); it != ov.end(); ++it )
      if ( !member(refs, *it) ) addReferences(*it, refs);
  }
}

void BaseRepository::
addInterfaces(const ClassDescriptionBase & db,
	      InterfaceMap & interfaceMap, bool all) {
  for ( ClassDescriptionBase::DescriptionVector::const_iterator it =
	  db.descriptions().begin(); it != db.descriptions().end(); ++it )
    if ( *it ) addInterfaces(**it, interfaceMap, all);
  TypeInterfaceMap::const_iterator cit = interfaces().find(&db);
  if ( cit == interfaces().end() ) return;
  for ( InterfaceSet::const_iterator iit = (cit->second).begin();
	iit != (cit->second).end(); ++iit ) {
    string n = (**iit).name();
    while ( all && member(interfaceMap, n) ) n = "+" + n;
    interfaceMap[n] = *iit;
  }
}

InterfaceMap BaseRepository::getInterfaces(const type_info & ti, bool all) {
  InterfaceMap interfaceMap;
  const ClassDescriptionBase * db = DescriptionList::find(ti);
  if ( !db ) return interfaceMap;
  addInterfaces(*db, interfaceMap, all);
  return interfaceMap;
}

void BaseRepository::
rebind(InterfacedBase & i, const TranslationMap & trans,
       const IVector & defaults) {
  InterfaceMap interfaceMap = getInterfaces(typeid(i), true);
  for ( InterfaceMap::iterator iit = interfaceMap.begin();
	iit != interfaceMap.end(); ++iit )
    iit->second->rebind(i, trans, defaults);
  i.rebind(trans);
}

void BaseRepository::update() {
  for_each(allObjects(), std::mem_fn(&InterfacedBase::update));
  clearAll(allObjects());
}

template <typename Set1, typename Set2>
bool overlap(const Set1 & s1, const Set2 & s2) {
  typename Set1::const_iterator i1 = s1.begin();
  typename Set2::const_iterator i2 = s2.begin();
  while ( i1 != s1.end() && i2 != s2.end() ) {
    if ( *i1 == *i2 ) return true;
    if ( *i1 < *i2 ) {
      i1 = s1.lower_bound(*i2);
      if ( *i1 == *i2 ) return true;
      ++i1;
    } else {
      i2 = s2.lower_bound(*i1);
      if ( *i1 == *i2 ) return true;
      ++i2;
    }
  }
  return false;
}

void BaseRepository::remove(tIBPtr ip) {
  ObjectMap::iterator it = objects().find(ip->fullName());
  if ( it == objects().end() || ip != it->second ) return;
  objects().erase(it);
  allObjects().erase(ip);
}

string BaseRepository::remove(const ObjectSet & rmset) {
  ObjectSet refset;
  for ( ObjectSet::const_iterator oi = rmset.begin();
	oi != rmset.end(); ++oi ) {
    IVector ov = GetObjectsReferringTo(*oi);
    refset.insert(ov.begin(), ov.end());
  }
  for ( ObjectSet::iterator oi = rmset.begin(); oi != rmset.end(); ++oi )
    refset.erase(*oi);
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
   
void BaseRepository::rename(tIBPtr ip, string newName) {
  ObjectSet::iterator it = allObjects().find(ip);
  if ( it == allObjects().end() ) {
    Register(ip, newName);
    return;
  }
  ObjectMap::iterator mit = objects().find(ip->fullName());
  if ( mit == objects().end() || mit->second != ip )
    throw RepoNameException(ip->fullName());
  
  objects().erase(mit);
  ip->name(newName);
  while ( member(objects(), ip->fullName()) ) ip->name(ip->fullName() + "#");
  objects()[ip->fullName()] = ip;
}

const InterfaceBase * BaseRepository::FindInterface(IBPtr ip, string name) {
  const auto & tmp=*ip;
  InterfaceMap imap = getInterfaces(typeid(tmp), false);
  InterfaceMap::iterator it = imap.find(name);
  return it == imap.end()? 0: it->second;
}

const ClassDocumentationBase * BaseRepository::getDocumentation(tcIBPtr ip) {
  const auto & tmp=*ip;
  TypeDocumentationMap::const_iterator cdoc =
    documentations().find(DescriptionList::find(typeid(tmp)));
  return cdoc != documentations().end()? cdoc->second: 0;
}

string BaseRepository::getModelDescription(tcIBPtr ip) {
  const ClassDocumentationBase *  cd = getDocumentation(ip);
  return cd? cd->modelDescription(): string("");
}

string BaseRepository::getModelReferences(tcIBPtr ip) {
  const ClassDocumentationBase *  cd = getDocumentation(ip);
  return cd? cd->modelReferences(): string("");
}

IBPtr BaseRepository::TraceObject(string path) {
  DirectoryAppend(path);
  string::size_type colon = path.find(':');
  IBPtr ip = GetPointer(path.substr(0, colon));
  if ( !ip ) {
    // Do special check if this is a decay mode.
    string name = path.substr(0, colon);
    string::size_type slash = name.rfind('/');
    if ( slash != string::npos ) name = name.substr(slash + 1);
    if ( name.find("->") != string::npos && name[name.length() - 1] == ';' ) {
      vector<DMPtr> save;
      DMPtr dm = DecayMode::constructDecayMode(name, &save);
      if ( dm )	ip = dynamic_ptr_cast<DMPtr>(GetPointer(path.substr(0, slash + 1)
							+ dm->tag()));
      if ( ip ) Throw<Exception>()
		  << "Warning: rewriting DecayMode name '"
		  << path.substr(0, colon).substr(slash + 1) << "' to '"
		  << ip->name() << Exception::warning;
    }
  }
  while ( colon != string::npos ) {
  if ( !ip ) throw RepositoryNotFound(path);
    path = path.substr(colon+1);
    colon = path.find(':');
    string::size_type bra = path.find('[');
    const InterfaceBase * ifb =
      FindInterface(ip, path.substr(0, min(colon, bra)));
    const ReferenceBase * rb = dynamic_cast<const ReferenceBase *>(ifb);
    if ( rb ) {
      ip = rb->get(*ip);
      continue;
    }
    const RefVectorBase * rvb = dynamic_cast<const RefVectorBase *>(ifb);
    if ( rvb ) {
      unsigned int place = 0;
      if ( bra < colon ) {
	string::size_type ket = path.find(']');
	place = atoi(path.substr(bra + 1,ket - bra - 1).c_str());
      }
      IVector iv = rvb->get(*ip);
      if ( place >= iv.size() ) throw RepositoryNotFound(path);
      ip = iv[place];
      continue;
    }
    const CommandBase * cb = dynamic_cast<const CommandBase *>(ifb);
    if ( cb ) {
      string::size_type ket = path.find(']');
      string newobj = cb->cmd(*ip, path.substr(bra + 1,ket - bra - 1));
      ip = GetPointer(newobj);
      continue;
    } 
    throw RepositoryNotFound(path);
  }
  if ( !ip ) throw RepositoryNotFound(path);
  return ip;
}

IBPtr BaseRepository::getObjectFromNoun(string noun) {
  string::size_type colon = noun.rfind(':');
  return TraceObject(noun.substr(0, colon));
}

string BaseRepository::getInterfaceFromNoun(string noun) {
  string::size_type colon = noun.rfind(':');
  string interface = noun.substr(colon+1);
  string::size_type bra = interface.find('[');
  if ( bra != string::npos ) return interface.substr(0, bra);
  else return interface;  
}

string BaseRepository::getPosArgFromNoun(string noun) {
  string::size_type colon = noun.rfind(':');
  string interface = noun.substr(colon+1);
  string::size_type bra = interface.find('[');
  if ( bra != string::npos ) {
    string::size_type ket = interface.find(']');
    return interface.substr(bra + 1,ket - bra - 1);
  }
  return "";
}

string BaseRepository::
GetInterfacedBaseClasses(const ClassDescriptionBase * cdb) {
  if ( !cdb || cdb->name() == "ThePEG::Interfaced" ||
       cdb->name() == "ThePEG::InterfacedBase" ) return "";
  string ret = cdb->name() + "\n";
  for ( int i = 0, N = cdb->descriptions().size(); i < N; ++i )
    ret += GetInterfacedBaseClasses(cdb->descriptions()[i]);
  return ret;
}

struct InterfaceOrder {
  bool operator()(const InterfaceBase * x, const InterfaceBase * y) const {
    return x->rank() > y->rank() ||
      ( x->rank() == y->rank() && x->name() < y->name() );
  }
};

void BaseRepository::readSetup(tIBPtr ip, istream & is) {
  ip->setup(is);
}

string BaseRepository::exec(string command, ostream &) {
  string verb = StringUtils::car(command);
  command = StringUtils::cdr(command);
  if ( verb.empty() || verb[0] == '#' ) return "";
  try {
    if ( verb == "DISABLEREADONLY" ) {
      InterfaceBase::NoReadOnly = true;
      return "";
    }
    if ( verb == "ENABLEREADONLY" ) {
      InterfaceBase::NoReadOnly = false;
      return "";
    }
    if ( verb == "cd" || verb == "pushd" || verb == "mkdir") {
      string dir = StringUtils::car(command);
      if ( verb == "cd" )
	ChangeDirectory(dir);
      else if ( verb == "pushd" )
	PushDirectory(dir);
      else
	CreateDirectory(dir);
      return "";
    }
    if ( verb == "popd" ) {
      PopDirectory();
      return "";
    }
    if ( verb == "pwd" ) return directoryStack().back();
    if ( verb == "dirs" ) {
      string ret;
      for ( StringVector::reverse_iterator it = directoryStack().rbegin();
	    it != directoryStack().rend(); ++it ) ret += *it;
      return ret;
    }
    if ( verb == "cp" || verb == "mv" ) {
      string oldname = StringUtils::car(command);
      DirectoryAppend(oldname);
      IBPtr obj = GetPointer(oldname);
      if ( !obj ) return "Error: No object named '" + oldname + "' available.";
      command = StringUtils::cdr(command);
      string newname = StringUtils::car(command);
      DirectoryAppend(newname);
      if ( newname[newname.size() - 1] == '/' ) newname += obj->name();
      if ( verb == "cp" ) obj = obj->fullclone();
      rename(obj, newname);
      return "";
    }
    if ( verb == "check" ) {
      string name = StringUtils::car(command);
      if ( directories().find(name) != directories().end() ) return name;
      if ( objects().find(name) != objects().end() ) return name;
      return "Not found";
    }
    if ( verb == "ls" ) {
      string className;
      string dir = StringUtils::car(command);
      if ( dir.size() ) {
	PushDirectory(dir);
	command = StringUtils::cdr(command);
	className = StringUtils::car(command);
      }
      string ret;
      string thisdir = directoryStack().back();
      for ( DirectorySet::iterator it = directories().begin();
	    it != directories().end(); ++it ) {
	string d = *it;
	if ( d.size() <= thisdir.size() ) continue;
	string d0 = d.substr(0, thisdir.size());
	string d1 = d.substr(thisdir.size());
	if ( d0 == thisdir && d1.find('/') == d1.size() - 1 ) {
	  if ( className.size() && SearchDirectory(d, className).empty() )
	    continue;
	  ret += (dir.size()? d: d1) + "\n";
	}
      }
      for ( ObjectMap::iterator it = objects().begin();
	    it != objects().end(); ++it ) {
	if ( className.size() ) {
	  const ClassDescriptionBase * cdb = DescriptionList::find(className);
    const auto & tmp=*(it->second);
	  if ( cdb &&
	       !DescriptionList::find(typeid(tmp))->isA(*cdb) )
	    continue;
	}
	if ( thisdir + it->second->name() == it->first )
	  ret += (dir.size()? it->first: it->second->name()) + '\n';
      }
      if ( dir.size() ) PopDirectory();
      return ret;
    }
    if ( verb == "library" ) {
      string library = StringUtils::car(command);
      if ( library.empty() ) return "Error: No library specified.";
      if ( !DynamicLoader::load(library) )
	return "Error: Could not load library " + library +
	  "\n - " + DynamicLoader::lastErrorMessage;
      return "";
    }

    if ( verb == "globallibrary" ) {
      string library = StringUtils::car(command);
      if ( library.empty() ) return "Error: No library specified.";
      if ( !DynamicLoader::load(library) )
	return "Error: Could not load library " + library +
	  "\n - " + DynamicLoader::lastErrorMessage;
      globalLibraries().push_back(library);
      return "";
    }

    if ( verb == "rmgloballibrary" ) {
      string library = StringUtils::car(command);
      if ( library.empty() ) return "Error: No library specified.";
      vector<string>::iterator it;
      while ( (it = find(globalLibraries(), library)) != globalLibraries().end() )
	globalLibraries().erase(it);
      return "";
    }

    if ( verb == "appendpath" ) {
      string path = StringUtils::car(command);
      if ( !path.empty() ) DynamicLoader::appendPath(path);
      return "";
    }

    if ( verb == "lspaths" ) {
      string paths;
      for ( int i = 0, N = DynamicLoader::allPaths().size(); i < N; ++i )
	paths += DynamicLoader::allPaths()[i] + "\n";
      return paths;
    }

    if ( verb == "prependpath" ) {
      string path = StringUtils::car(command);
      if ( !path.empty() ) DynamicLoader::prependPath(path);
      return "";
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
      if ( !db )  {
	string msg = "Error: " + className + ": No such class found.";
	if ( !DynamicLoader::lastErrorMessage.empty() )
	  msg += "\nerror message from dynamic loader:\n" +
	    DynamicLoader::lastErrorMessage;
	return msg;
      }
      IBPtr obj = dynamic_ptr_cast<IBPtr>(db->create());
      if ( !obj ) return "Error: Could not create object of class "+className;
      if ( name.empty() ) return "Error: No name specified.";
      Register(obj, name);
      return "";
    }
    if ( verb == "setup" ) {
      string name = StringUtils::car(command);
      DirectoryAppend(name);
      IBPtr obj = GetPointer(name);
      if ( !obj ) return "Error: Could not find object named " + name;
      istringstream is(StringUtils::cdr(command));
      readSetup(obj, is);
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
    if ( verb == "rcp" ) {
      string name = StringUtils::car(command);
      DirectoryAppend(name);
      string newName = StringUtils::car(StringUtils::cdr(command));
      if ( newName.empty() )
	return "Error: No destination directory specified.";
      DirectoryAppend(newName);
      CreateDirectory(newName);
      if ( newName[newName.size() - 1]  != '/' ) newName += '/';
      IBPtr obj = GetPointer(name);
      if ( name[name.size() - 1]  != '/' ) name += '/';
      IVector ov = SearchDirectory(name);
      ov.push_back(obj);
      if ( ov.empty() ) return "Error: No such object or directory.";
      ObjectSet toclone;
      for ( IVector::iterator i = ov.begin(); i != ov.end(); ++i ) {
	toclone.insert(*i);
	addReferences(*i, toclone);
      }
      for ( ObjectSet::iterator i = toclone.begin(); i != toclone.end(); ++i )
	Register((**i).clone(), newName + (**i).name());
      return "";
    }
    if ( verb == "rebind" ) {
      // For all objects in the repository, replace any references to
      // the first object given with references to the second
      // object. The two objects will change names
     IBPtr ip1 = TraceObject(StringUtils::car(command));
     string newname = StringUtils::car(StringUtils::cdr(command));
     DirectoryAppend(newname);
     IBPtr ip2 = GetPointer(newname);
     if ( !ip2 ) {
       ip2 = ip1->fullclone();
       rename(ip2, newname);
     }
     TranslationMap trans;
     trans[ip1] = ip2;
     IVector objs = GetObjectsReferringTo(ip1);
     for ( int i = 0, N = objs.size(); i < N; ++i ) rebind(*objs[i], trans, IVector());
    }
    if ( verb == "doxygendump" ) {
      string spacename = StringUtils::car(command);
      command = StringUtils::cdr(command);
      string filename = StringUtils::car(command);
      ofstream os(filename.c_str());
      for ( TypeDocumentationMap::const_iterator it = documentations().begin();
	    it != documentations().end(); ++it ) {
	const ClassDescriptionBase & db = *(it->first);
	string classname = db.name();
	if ( classname.substr(0, spacename.length()) != spacename ) continue;
	string briefname = classname.substr(spacename.length());
	os << "/** \\page " << briefname << "Interfaces "
	   << "Interfaces defined for the " << classname << " class.\n\n"
	   << "\\par Brief class description:\n";
	string doc = it->second->documentation();
	if ( doc.substr(0,25) == "There is no documentation" )
	  os << "See " << classname << "\n\n";
	else
	  os << doc << "<br>See also " << classname << "\n\n";
	TypeInterfaceMap::const_iterator isit = interfaces().find(it->first);
	if ( isit == interfaces().end() || isit->second.empty() ) {
	  os << "There are no interfaces declared for this class.\n\n";
	} else {
	  const InterfaceSet & ints = isit->second;
	  for ( InterfaceSet::const_iterator iit = ints.begin();
		iit != ints.end(); ++iit )
	    (**iit).doxygenDescription(os);
	}
	string baserefs = "";
	int nbases = 0;
	for ( int ib = 0, N = db.descriptions().size(); ib < N; ++ib ) {
	  if ( documentations().find(db.descriptions()[ib]) ==
	       documentations().end() ) continue;
	  const ClassDescriptionBase & bdb = *db.descriptions()[ib];
	  if ( nbases ) baserefs += " and ";
	  string briefname = bdb.name().substr(bdb.name().rfind("::") + 2);
	  baserefs += "\\ref " + briefname +
	    "Interfaces \"" + bdb.name() + "\"";
	  ++nbases;
	}
	if ( nbases == 1 )
	  os << "<hr>There may be interfaces inherited from the "
	     << baserefs << " class.";
	else if ( nbases > 1 ) 
	  os << "<hr>There may be interfaces inherited from the "
	     << "following classes: " << baserefs << ".";
	os << "\n\n*/\n\n";
      }
      return "";
    }
    if ( verb == "mset" || verb == "msetdef" || verb == "minsert" ||
	 verb == "mdo" || verb == "mget" || verb == "mdef" || verb == "mmin" ||
	 verb == "mmax" || verb == "merase" || verb == "msend"  ) {
      if ( verb == "msend" ) verb = "mdo";
      string dir = StringUtils::car(command);
      command = StringUtils::cdr(command);
      string className = StringUtils::car(command);
      command = StringUtils::cdr(command);
      string interface = StringUtils::car(command);
      string arguments = StringUtils::cdr(command);
      string::size_type bra = interface.find('[');
      if ( bra != string::npos ) {
	string::size_type ket = interface.find(']');
	arguments = interface.substr(bra + 1,ket - bra - 1) + " " + arguments;
	interface = interface.substr(0, bra);
      }
      IVector ov = SearchDirectory(dir, className);
      if ( ov.empty() ) return "Error: no matching objects found.";
      string ret;
      verb = verb.substr(1);
      for ( IVector::size_type i = 0; i < ov.size(); ++i ) {
	const InterfaceBase * ifb = FindInterface(ov[i], interface);
	if ( !ifb ) continue;
	string mess = ifb->exec(*ov[i], verb, arguments);
	if ( !mess.empty() ) ret += ov[i]->fullName() + ": " + mess + "\n";
      }
      return ret.substr(0, ret.size() - 1);
    }
    if ( verb == "set" || verb == "setdef" || verb == "insert" ||
	 verb == "do" || verb == "get" || verb == "def" || verb == "min" ||
	 verb == "max" || verb == "describe" || verb == "fulldescribe" ||
	 verb == "erase" || verb == "clear" || verb == "send" || verb == "newdef" ) {
      if ( verb == "send" ) verb = "do";
      if ( verb == "newdef" && !InterfaceBase::NoReadOnly )
	return "Error: The default value of an interface is a read-only "
	  "entity. Use the command 'DISABLEREADONLY' to override.";
      string noun = StringUtils::car(command);
      string arguments = getPosArgFromNoun(noun) + " "
	+ StringUtils::cdr(command);
      IBPtr ip = getObjectFromNoun(noun);
      const InterfaceBase * ifb = FindInterface(ip, getInterfaceFromNoun(noun));
      if ( !ifb && verb != "describe" && verb != "fulldescribe" ) {
	string ret = "Error: The interface '" + noun + "' was not found.\n";
	ret += "Valid interfaces:\n";
  const auto & tmp=*ip;
	InterfaceMap imap = getInterfaces(typeid(tmp));
	for ( InterfaceMap::iterator it = imap.begin(); it != imap.end(); ++it )
	  ret += "* " + it->second->name() + "\n";
	return ret;
      }
      if ( verb == "describe" ) {
	if ( ifb ) return ifb->description();
  const auto & tmp=*ip;
	const ClassDescriptionBase * cd = DescriptionList::find(typeid(tmp));
	string ret = "Object '" + ip->name() + "' of class '" +
	  cd->name() + "':\n";
	TypeDocumentationMap::const_iterator cdoc = documentations().find(cd);
	if ( cdoc != documentations().end() )
	  ret += cdoc->second->documentation() + "\n";
	ret +="Interfaces:\n";
	InterfaceMap imap = getInterfaces(typeid(tmp));
	for ( InterfaceMap::iterator it = imap.begin(); it != imap.end(); ++it )
	  ret += "* " + it->second->name() + "\n";
	return ret;
      } else if ( verb == "fulldescribe" ) {
	if ( ifb ) return ifb->fullDescription(*ip);
	ostringstream ret;
  const auto & tmp=*ip;
	const ClassDescriptionBase * cd = DescriptionList::find(typeid(tmp));
	TypeDocumentationMap::const_iterator cdoc = documentations().find(cd);
	ret << ip->fullName() << endl << cd->name() << endl;
	if ( cdoc != documentations().end() )
	  ret << cdoc->second->documentation() << endl;
	ret << "Interfaces:" << endl;
	InterfaceMap imap = getInterfaces(typeid(tmp));
	typedef set<const InterfaceBase *, InterfaceOrder> InterfaceSet;
	InterfaceSet iset;
	for ( InterfaceMap::iterator it = imap.begin(); it != imap.end(); ++it )
	  iset.insert(it->second);
	double rank = 1.0;
	for ( InterfaceSet::iterator it = iset.begin();
	      it != iset.end(); ++it ) {
	  if ( rank >= 0.0 && (**it).rank() < 0.0 ) ret << "0" << endl;
	  rank = (**it).rank();
	  ret << (**it).type() << " " << (**it).name() << endl;
	}
	return ret.str();
      }
      else
	return ifb->exec(*ip, verb, arguments);
    }
    if ( verb == "baseclasses" ) {
      string className = StringUtils::car(command);
      const ClassDescriptionBase * cdb = 0;
      if ( className.size() ) {
	cdb = DescriptionList::find(className);
	if ( !cdb ) return "Error: no class '" + className + "' found.";
      }
      return GetInterfacedBaseClasses(cdb);
    }
    if ( verb == "describeclass" ) {
      string className = StringUtils::car(command);
      const ClassDescriptionBase * cdb = 0;
      if ( className.size() ) {
	cdb = DescriptionList::find(className);
	if ( !cdb ) return "Error: no class '" + className + "' found.";
      }
      TypeDocumentationMap::const_iterator cdoc = documentations().find(cdb);
      if ( cdoc != documentations().end() )
	return cdoc->second->documentation() + "\n";
      else
	return "";
    }
    if ( verb == "lsclass" ) {
      string className = StringUtils::car(command);
      const ClassDescriptionBase * cdb = 0;
      if ( className.size() ) {
	cdb = DescriptionList::find(className);
	if ( !cdb ) return "Error: no class '" + className + "' found.";
      }
      vector<const ClassDescriptionBase *> classes;
      if ( cdb && !cdb->abstract() ) classes.push_back(cdb);
      for ( DescriptionList::DescriptionMap::const_iterator
	      it = DescriptionList::all().begin();
	    it != DescriptionList::all().end(); ++it ) {
	if ( it->second == cdb || it->second->abstract() ) continue;
	if ( cdb && !it->second->isA(*cdb) ) continue;
	classes.push_back(it->second);
      }
      if ( classes.empty() )  return "Error: no classes found.";
      string ret;
      for ( int i = 0, N = classes.size(); i < N; ++i ) 
	ret += classes[i]->name() + "\n";
      return ret;
    }
    
  }
  catch (const Exception & e) {
    e.handle();
    return "Error: " + e.message();
  }
  return "Error: Unrecognized command '" + verb + "'.";
}

BadClassClone::BadClassClone(const InterfacedBase & o) {
  theMessage << "Could not clone the object '" << o.name()
	     << "' of class '" << TypeInfo::name(o)
	     << "' because the class does not"
	     << " implement a working 'clone' method.";
  severity(abortnow);
}

BadClone::BadClone(const InterfacedBase & o) {
  theMessage << "Could not clone the object '" << o.name()
	     << "' of class '" << TypeInfo::name(o)
	     << "' because the clone method threw an unknown exception.";
  severity(abortnow);
}

RepoNameException::RepoNameException(string name) {
  theMessage << "The object '" << name << "' is present in the Repository but "
	     << "under a different name. This means that the name of the "
	     << "object has been illegally changed outside of the Repository.";
  severity(abortnow);
}

RepoNameExistsException::RepoNameExistsException(string name) {
  theMessage
    << "The object '" << name
    << "' was not created as another object with that name already exists.";
  severity(warning);
}

RepositoryNoDirectory::RepositoryNoDirectory(string name) {
  theMessage << "The directory '" << name << "' does not exist.";
  severity(warning);
}

RepositoryNotFound::RepositoryNotFound(string name) {
  theMessage << "There was no object named '" << name << "' in the repository.";
  severity(warning);
}

RepositoryClassMisMatch::
RepositoryClassMisMatch(const InterfacedBase & o, string name) {
  theMessage << "The requested object '" << o.fullName() << "' was not of the "
	     << "specified type (" << name << ").";
  severity(warning);
}

