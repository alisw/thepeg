// -*- C++ -*-
//
// Tree.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_Tree_H
#define LWH_Tree_H
//
// This is the declaration of the Tree class.
//

#include "AITree.h"
#include "ManagedObject.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>

namespace LWH {

using namespace AIDA;

/**
 * The Tree class is a simple implementation of the AIDA::ITree
 * interface.
 */
class Tree: public ITree {

public:

  /** The AnalysisFactory is a friend. */
  friend class AnalysisFactory;

  /** A path is a vector of directory names. */
  typedef std::vector<std::string> Path;

  /** A set of paths */
  typedef std::set<Path> PathSet;

  /** Map of paths to objects. */
  typedef std::map<std::string, IManagedObject *> ObjMap;

public:

  /**
   * The standard constructor.
   */
  Tree(std::string storename, bool xml = true)
    : name(storename), flat(!xml), cwd("/"), overwrite(true) {
    dirs.insert(Path());
  }

  /**
   * The default constructor.
   */
  Tree(): name(""), flat(false), cwd("/") {
    dirs.insert(Path());
  }

  /**
   * The copy constructor.
   */
  Tree(const Tree & dt)
    : ITree(dt), name(dt.name), flat(dt.flat), dirs(dt.dirs),
      objs(dt.objs), cwd(dt.cwd), overwrite(true) {}

  /// Destructor.
  virtual ~Tree() {
    for ( ObjMap::iterator it = objs.begin(); it != objs.end(); ++it )
      delete it->second;
  }

  /**
   * Get the name of the store.
   * @return The store's name.
   */
  std::string storeName() const {
    return name;
  }

  /**
   * Get the IManagedObject at a given path in the ITree. The path can either be
   * absolute or relative to the current working directory.
   * @param path The path.
   * @return     The corresponding IManagedObject.
   */
  IManagedObject * find(const std::string & path) {
    ObjMap::const_iterator it = objs.find(path);
    return it == objs.end()? (IManagedObject *)0: it->second;
  }

  /**
   * LWH cannot get a mounted ITree at a given path in the current ITree.
   * @return     0 always.
   */
  ITree * findTree(const std::string &) {
    return 0;
  }

  /**
   * Change to a given directory.
   * @param dir The absolute or relative path of the directory we are
   * changing to.
   * @return false If the path does not exist.
   */
  bool cd(const std::string & dir) {
    PathSet::iterator it = dirs.find(purgepath(str2pth(fullpath(sts(dir)))));
    if ( it == dirs.end() ) return false;
    cwd = pth2str(*it);
    return true;
  }

  /**
   * Insert the ManagedObject \a o in the tree with the path \a str.
   */
  bool insert(std::string str, IManagedObject * o) {
    Path path = purgepath(str2pth(fullpath(str)));
    if ( dirs.find(path) == dirs.end() ) {
      std::string fullname = pth2str(path);
      path.pop_back();
      if ( dirs.find(path) != dirs.end() ) {
	ObjMap::iterator old = objs.find(fullname);
	if ( old == objs.end() || overwrite ) {
	  if ( old != objs.end() ) {
	    delete old->second;
	    objs.erase(old);
	  }
	  objs[fullname] = o;
	  return true;
	}
      }
    }
    return false;
  }

  /**
   * Get the path of the current working directory.
   * @return The path of the current working directory.
   */
  std::string pwd() const {
    return cwd;
  }

  /** 
   * Not implemented in LWH.
   * @return false always.
   *
   */
  bool ls(const std::string & = ".", bool = false,
	  std::ostream & = std::cout) const {
    return false;
  }

  /**
   * Not implemented in LWH.
   */
  std::vector<std::string> listObjectNames(const std::string & = ".",
					   bool = false) const {
    return std::vector<std::string>();
  }

  /**
   * Not implemented in LWH.
   */
  std::vector<std::string> listObjectTypes(const std::string & = ".",
					   bool = false) const {
    return std::vector<std::string>();
  }

  /**
   * Create a new directory. Given a path only the last directory
   * in it is created if all the intermediate subdirectories already exist.
   * @param dir The absolute or relative path of the new directory.
   * @return false If a subdirectory within the path does
   * not exist or it is not a directory. Also if the directory already exists.
   */   
  bool mkdir(const std::string & dir) {
    Path p = purgepath(str2pth(fullpath(sts(dir))));
    Path base = p;
    base.pop_back();
    if ( dirs.find(base) == dirs.end() ) return false;
    dirs.insert(p);
    return true;
  }

  /**
   * Create a directory recursively. Given a path the last directory
   * and all the intermediate non-existing subdirectories are created.
   * @param dir The absolute or relative path of the new directory.
   * @return false If an intermediate subdirectory
   *             is not a directory, or if the directory already exists.
   */
  bool mkdirs(const std::string & dir) {
    return mkdirs(purgepath(str2pth(fullpath(sts(dir)))));
  }

  /**
   * Create a directory recursively. Given a Path the last directory
   * and all the intermediate non-existing subdirectories are created.
   * @param p The full Path of the new directory.
   * @return false If an intermediate subdirectory
   *             is not a directory, or if the directory already exists.
   */
  bool mkdirs(Path p) {
    if ( dirs.find(p) != dirs.end() ) return true;
    dirs.insert(p);
    p.pop_back();
    return mkdirs(p);
  }

  /**
   * Remove a directory and all the contents underneeth.
   * @param dir The absolute or relative path of the directory to be removed.
   * @return false If path does not exist or if it is not
   *             a directory or if the directory is not empty.
   */
  bool rmdir(const std::string & dir) {
    Path path = purgepath(str2pth(fullpath(sts(dir))));
    if ( dirs.find(path) == dirs.end() ) return false;
    for ( ObjMap::const_iterator it = objs.begin(); it != objs.end(); ++it )
      if ( it->first.substr(0, dir.length()) == dir ) return false;
    dirs.erase(path);
    return true;
  }

  /**
   * Remove and delete an IManagedObject by specifying its path.
   * @param path The absolute or relative path of the IManagedObject to be
   * removed.
   * @return false If path does not exist.
   */
  bool rm(const std::string & path) {
    ObjMap::iterator it = objs.find(fullpath(path));
    if ( it == objs.end() ) return false;
    delete it->second;
    objs.erase(it);
    return true;
  }

  /**
   * Get the full path of an IManagedObject.
   * @param o The IManagedObject whose path is to be returned.
   * @return  The object's absolute path.
   *          If the object does not exist, an empty string is returned.
   */
  std::string findPath(const IManagedObject & o) const {
    for ( ObjMap::const_iterator it = objs.begin(); it != objs.end(); ++it )
      if ( it->second == &o ) return it->first;
    return "";
  }

  /**
   * Move an IManagedObject or a directory from one directory to another.
   * @param oldp The path of the IManagedObject [not direcoty] to be moved.
   * @param newp The path of the diretory in which the object has to be
   * moved to.
   * @return false If either path does not exist.
   */
  bool mv(const std::string & oldp, const std::string & newp) {
    Path newpath = purgepath(str2pth(fullpath(sts(newp))));
    std::string foldp = fullpath(oldp);
    Path oldpath = purgepath(str2pth(foldp));
    ObjMap::iterator it = objs.find(foldp);
    if ( it == objs.end() ) return false;
    if ( dirs.find(newpath) != dirs.end() ) return false;
    newpath.push_back(oldpath.back());
    if ( !insert(pth2str(newpath), it->second) ) return false;
    objs.erase(foldp);
    return true;
  }

  /**
   * Print all histograms to the current filename.
   * @return false if something went wrong.
   */
  bool commit() {
    std::ofstream of(name.c_str());
    if ( !of ) return false;
    if ( !flat ) of
      << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE aida SYSTEM "
      << "\"http://aida.freehep.org/schemas/3.0/aida.dtd\">\n"
      << "<aida version=\"3.0\">\n"
      << "<implementation version=\"1.0\" package=\"FreeHEP\"/>" << std::endl;
    for ( ObjMap::const_iterator it = objs.begin(); it != objs.end(); ++it ) {
      ManagedObject * o = dynamic_cast<ManagedObject *>(it->second);
      if ( !o ) continue;
      std::string path = it->first.substr(0, it->first.rfind('/'));
      std::string name = it->first.substr(it->first.rfind('/') + 1);
      if ( flat )
	o->writeFLAT(of, path, name);
      else
	o->writeXML(of, path, name);
    }
    if ( !flat ) of << "</aida>" << std::endl;
    return of.good();
  }

  /**
   * Not implemented in LWH.
   */
  void setOverwrite(bool o = true) {
    overwrite = o;
  }

  /**
   * Not implemented in LWH.
   * @return false always.
   */
  bool cp(const std::string &, const std::string &, bool = false) {
    return false;
  }

  /**
   * Not implemented in LWH.
   * @return false always.
   */
  bool symlink(const std::string &, const std::string &) {
    return false;
  }

  /**
   * Not implemented in LWH.
   * @return false always.
   */
  bool mount(const std::string &, ITree &, const std::string &) {
    return false;
  }

  /**
   * Not implemented in LWH.
   * @return false always.
   */
  bool unmount(const std::string &) {
    return false;
  }

  /**
   * Calls commit().
   */
  bool close() {
    return commit();
  }

  /**
   * Not implemented in LWH.
   * @return null pointer always.
   */ 
  void * cast(const std::string &) const {
    return 0;
  }

protected:

  /** Strip trailing slash. */
  std::string sts(std::string s) const {
    if ( s[s.length() - 1] == '/' ) s = s.substr(0, s.length() - 1);
    if ( s[s.length() - 1] == '/' ) return "";
    return s;
  }

  /** Strip trailing name */
  std::string stn(std::string s) const {
    std::string::size_type slash = s.rfind('/');
    return s.substr(0, slash);
  }

  /** Get proper full path from possibly relative path. */
  std::string fullpath(std::string d) const {
    if ( d[0] != '/' ) d = cwd + "/" + d;
    return pth2str(purgepath(str2pth(d)));
  }

  /** Convert a string containing a path to a Path object. */
  Path str2pth(std::string s) const {
    Path pth;
    std::string::size_type i = s.find_first_not_of("/");
    while ( i != std::string::npos ) {
      s = s.substr(i);
      i = s.find_first_of("/");
      pth.push_back(s.substr(0, i));
      if ( i == std::string::npos ) return pth;
      s = s.substr(i);
      i = s.find_first_not_of("/");
    }
    return pth;
  }

  /** Convert a Path object to a corresponding string. */
  std::string pth2str(const Path & pth) const {
    std::string str;
    for ( int i = 0, N = pth.size(); i < N; ++i ) str += "/" + pth[i];
    return str;
  }

  /** Remove '..' and '.' components of the given Path object. */
  Path purgepath(const Path & pth) const {
    Path p;
    for ( int i = 0, N = pth.size(); i < N; ++i ) {
      if ( pth[i] == ".." ) p.pop_back();
      else if ( pth[i] != "." ) p.push_back(pth[i]);
    }
    return p;
  }

private:

  /** The filename to print histograms to. */
  std::string name;

  /** If true write histograms in FLAT format, otherwise in XML. */
  bool flat;

  /** The set of defined directories. */
  PathSet dirs;

  /** The set of defined objects. */
  ObjMap objs;

  /** The current working directory. */
  std::string cwd;

  /** Overwrite strategy. */
  bool overwrite;

};

}

#endif /* LWH_Tree_H */
