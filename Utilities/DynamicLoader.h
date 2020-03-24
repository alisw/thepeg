// -*- C++ -*-
//
// DynamicLoader.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_DynamicLoader_H
#define ThePEG_DynamicLoader_H
// This is the declaration of the DynamicLoader class.

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * <code>DynamicLoader</code> is the general interface to the dynamic
 * loader functions of the underlying operating system. Currently it
 * only works on Linux.
 *
 * @see ClassTraits
 * @see ClassDescription
 * @see DescriptionList
 * @see PersistentIStream
 */
class DynamicLoader {

public:

  /**
   * The actual load command used on the current platform.
   */
  static bool loadcmd(string);

  /**
   * Try to load the file given as argument. If the filename does not
   * begin with a '/', try to prepend the paths one at the time until
   * success.  If all fail try without prepending a path.
   * @return true if the loading succeeded, false otherwise.
   */
  static bool load(string file);

  /**
   * Add a path to the bottom of the list of directories to seach for
   * dynaically linkable libraries.
   */
  static void appendPath(string);

  /**
   * Add a path to the top of the list of directories to seach for
   * dynaically linkable libraries.
   */
  static void prependPath(string);

  /**
   * Return the last error message issued from the platforms loader.
   */
  static string lastErrorMessage;

  /**
   * Insert the name of the given library with correct version numbers
   * appended, in the corresponding map.
   */
  static void dlname(string);

  /**
   * Given a list of generic library names, return the same list with
   * appended version numbers where available.
   */
  static string dlnameversion(string libs);

  /**
   * Return the full list of directories to seach for dynaically
   * linkable libraries.
   */
  static const vector<string> & allPaths();

  /**
   * Return the list of appended directories to seach for dynaically
   * linkable libraries.
   */
  static const vector<string> & appendedPaths();

  /**
   * Return the list of prepended directories to seach for dynaically
   * linkable libraries.
   */
  static const vector<string> & prependedPaths();

private:

  /**
   * The list of directories to seach for dynaically linkable
   * libraries.
   */
  static vector<string> paths;

  /**
   * The list of prepended directories to seach for dynaically linkable
   * libraries.
   */
  static vector<string> prepaths;

  /**
   * The list of appended directories to seach for dynaically linkable
   * libraries.
   */
  static vector<string> apppaths;

  /**
   * Used to initialize the paths vector from the ThePEG_PATH
   * environment.
   */
  static vector<string> defaultPaths();

  /**
   * Map of names of dynamic libraries with correct version numbers
   * indexed by their generic names.
   */
  static map<string,string> versionMap;

};

}

#endif /* ThePEG_DynamicLoader_H */
