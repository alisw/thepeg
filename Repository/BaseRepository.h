// -*- C++ -*-
//
// BaseRepository.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_BaseRepository_H
#define ThePEG_BaseRepository_H
// This is the declaration of the BaseRepository class.

#include "ThePEG/Config/ThePEG.h"
#include "BaseRepository.xh"
#include "ThePEG/Interface/InterfaceBase.fh"
#include "ThePEG/Interface/ClassDocumentation.fh"
#include "ThePEG/Interface/InterfacedBase.h"
#include "ThePEG/Utilities/ClassDescription.fh"

namespace ThePEG {

/**
 * BaseRepository is a purely static class which keeps a set of
 * InterfacedBase objects indexed by their name. The objects and their
 * names are divided up in a tree-like structure inspired by the Unix
 * file system.
 *
 * The InterfacedBase objects may be manipulated using InterfaceBase
 * objects. This may be done directly or via a simple command
 * interface using the exec() method.
 *
 * RepositoryBase is closely related to the Repository sub-class. The
 * division may seem unnecessary, but the idea is that BaseRepository
 * is a general repository for administrating and manipulating a set
 * of InterfacedBase objects, while the Repository adds on utilites
 * which are special to ThePEG where the objects are Interfaced (a
 * sub-class of InterfacedBase).
 *
 * @see Repository
 * @see InterfacedBase
 * @see InterfaceBase
 * @see Interfaced
 * 
 */
class BaseRepository {

public:

  /** A set of strings. */
  typedef StringSet DirectorySet;

  /** A vector of character strings. */
  typedef vector<string> StringVector;

  /** A set of pointers to InterfaceBase objects. */
  typedef set<const InterfaceBase *> InterfaceSet;

  /** A map of sets of IterfaceBase objects indexed by pointers to
      ClassDescriptionBase objects. */
  typedef map<const ClassDescriptionBase *, InterfaceSet> TypeInterfaceMap;

  /** A map of ClassDocumentationBase objects indexed by pointers to
      ClassDescriptionBase objects. */
  typedef map<const ClassDescriptionBase *, const ClassDocumentationBase *>
    TypeDocumentationMap;
 
public:

  /**
   * Interpret the command in \a cmd and return possible
   * messages. This is the main function for the command-line
   * interface. The syntax is described elsewhere. The ostream
   * argument is currently unused.
   */
  static string exec(string cmd, ostream &);

  /** @name Functions for adding and deleting objects and interfaces. */
  //@{
  /**
   * Register an interface. This is called automatically in the
   * InterfaceBase constructor and should never be called explicitly.
   */
  static void Register(const InterfaceBase &, const type_info &);

  /**
   * Register a class documentation. This is called automatically in
   * the ClassDocumentationBase constructor and should never be called
   * explicitly.
   */
  static void Register(const ClassDocumentationBase &, const type_info &);

  /**
   * Register a new object using the its current name. If the object
   * is already in the repository, nothing happens. If another object
   * already exists with the same name, the new object will have
   * <code>#</code>'s appended to its name to make it unique.
   */
  static void Register(IBPtr);

  /**
   * Register a new object giving it a new \a name. If the object is
   * already in the repository, nothing happens. If another object
   * already exists with the same name, the new object will have
   * <code>#</code>'s appended to its name to make it unique.
   */
  static void Register(IBPtr, string name);

  /**
   * Remove the given object from the repository. If the object was
   * not present nothing will happen.
   */
  static void remove(tIBPtr);

  /**
   * Remove objects. Remove the objects in \a rmset if there are no
   * other objects in the repository referring to them, otherwise
   * return an error message and the names of the objects refering to
   * them separated by new-line characters.
   */
  static string remove(const ObjectSet & rmset);

  /**
   * Rename a given \a object. Syntacticly the same as
   * <code>remove(object); Register(object, newName);</code>.
   */
  static void rename(tIBPtr object, string newName);
  //@}

  /** @name Access the directory stack. */
  //@{
  /**
   * Create a new directory with the given name. If the given name
   * starts with a <code>/</code> the name is assumed to be an absolute
   * path, otherwise it is assumed to be a path relative to the
   * current directory.
   */
  static void CreateDirectory(string);

  /**
   * Check if directory exixts. Check if the name given as argument
   * corresponds to an existing directory. If the argument string does
   * not end in a <code>/</code> it is assumed to be the name of an
   * object in a directory, and only the directory part of the name is
   * checked. If the given name starts with a <code>/</code> the name
   * is assumed to be an absolute path, otherwise it is assumed to be
   * a path relative to the current directory.
   *
   * @throws RepositoryNoDirectory if the correspinding directory is
   * non-existent.
   */
  static void CheckObjectDirectory(string);

  /**
   * Check if directory exixts. Check if the name given as argument
   * corresponds to an existing directory. If the given name starts
   * with a <code>/</code> the name is assumed to be an absolute path,
   * otherwise it is assumed to be a path relative to the current
   * directory.
   *
   * @throws RepositoryNoDirectory if the correspinding directory is
   * non-existent.
   */
  static void CheckDirectory(string);

  /**
   * Return the absolute path.  If the given name starts with a
   * <code>/</code> the name is assumed to be an absolute path already,
   * otherwise it is assumed to be a path relative to the current
   * directory, and the absolute path is constructed.
   */
  static void DirectoryAppend(string &);

  /**
   * Set the current directory to \a name. \a name can be aither a
   * relative or absolute path. The new directory replaces the
   * previous current directory on the directory stack.
   *
   * @throws RepositoryNoDirectory if the directory is non-existent.
   */
  static void ChangeDirectory(string name);

  /**
   * Set the current directory to \a name. \a name can be aither a
   * relative or absolute path. The new directory is pushed onto the
   * directory stack.
   *
   * @throws RepositoryNoDirectory if the directory is non-existent.
   */
  static void PushDirectory(string name);

  /**
   * Pop the directory stack. Leave the current directory and set the
   * directory which is on top of the popped directory stack.
   */
  static void PopDirectory();

  /**
   * A list of all globally loaded libraries.
   */
  static vector<string> & globalLibraries();

  //@}

  /** @name Information on where to read input files. */
  //@{
protected:

  /**
   * The stack of directories used by the "read" command.
   */
  static stack<string> & currentReadDirStack();

  /**
   * List of directories to search for files for the "read" command.
   */
  static vector<string> & readDirs();

public:

  /**
   * Access to list of directories to search for files for the "read" command.
   */
  static const vector<string> & listReadDirs();

  /**
   * Add a directory to readDirs().
   */
  static void prependReadDir(string);
  
    /**
   * Add a string vector with directories to readDirs().
   */
  static void prependReadDir(const std::vector<std::string>& dirs);

  /**
   * Add a directory to readDirs().
   */
  static void appendReadDir(string);
  
    /**
   * Add a string vector with directories to readDirs().
   */
  static void appendReadDir(const std::vector<std::string>& dirs);

  //@}

  /** @name Access objects in the repository. */
  //@{
  /**
   * Return a reference counted pointer to the given object. This
   * currently not needed when ThePEG is used with the
   * ThePEG::Pointer::RCPtr class of pointers.
   */
  template <typename T>
  static typename Ptr<T>::pointer GetPtr(const T &);

  /**
   * Return a pointer of the specified type to an object with the
   * given name. If such an object does not exist, GetPtr will return
   * a null pointer.
   */
  template <typename PtrType>
  static PtrType GetPtr(string);

  /**
   * Return a pointer of the specified type to an object with the
   * given name. If such an object does not exist an exception will be
   * thrown.
   * @throws RepositoryNotFound if the object was not found.
   * @throws RepositoryClassMisMatch if the object exists but is of
   * the wrong class.
   */
  template <typename PtrType>
  static PtrType GetObject(string);

  /**
   * Return a pointer to an object with the given name or null if no
   * such object exists.
   */
  static IBPtr GetPointer(string);

  /**
   * Return all objects in the directory \a name. Optionally only return
   * objects of class \a className or of a sub-class thereof.
   */
  static IVector SearchDirectory(string name, string className = "");

  /**
   * Find an object. If the \a name does not begin with '/', the
   * current directory is prepended. If the string is on the form
   * <code>object:interface</code> (or
   * <code>object:interface[i]</code>) and <code>interface</code>
   * corresponds to an Reference (or RefVector) interface, the
   * corresponding referenced object is returned. (also
   * <code>object:interface:interface</code> is allowed etc.)
   */
  static IBPtr TraceObject(string name);

  /**
   * Return a string containing the name of the given class
   * description and its base classes, one on each line.
   */
  static string GetInterfacedBaseClasses(const ClassDescriptionBase * cdb);

  /**
   * Get an object. Decompose a string of the form
   * <code>object:interface</code> or
   * <code>object:vector-interface[pos]</code>. Retrun a pointer to
   * the corresponding <code>object</code>.
   */
  static IBPtr getObjectFromNoun(string noun);
  //@}

  /** @name Access references between object in the repository. */
  //@{
  /**
   * Get referring objects. Return all object which refers to the
   * given object through a Reference of RefVector interface.
   */
  static IVector GetObjectsReferringTo(IBPtr);

  /**
   * Get direct references. Return all objects the given object refers
   * to directly through a Reference of RefVector interface.
   */
  static IVector DirectReferences(IBPtr);

  /**
   * Get all references. If \a obj contains references to other objects,
   * either through a Reference or RefVector interface or through the
   * virtual getReferences member function, add these to refs. Do the
   * same to the references recursively.
   */
  static void addReferences(tIBPtr obj, ObjectSet & refs);
  //@}

  /** @name Access the interfaces of the objects in the repository. */
  //@{
  /**
   * Get interfaces. Return the interfaces defined for the
   * InterfacedBase class with the given type_info, \a ti, mapped to
   * their name. If several interfaces with the same name exists only
   * the one which correspond to the most derived class will be given,
   * except if \a all is true in which case all interfaces are given
   * (prefixed by '+'s to become unique).
   */
  static InterfaceMap getInterfaces(const type_info & ti, bool all = true);

  /**
   * Return an interface with the given \a name to the given \a object.
   */
  static const InterfaceBase * FindInterface(IBPtr object, string name);

  /**
   * Get an interface name. Decompose a string of the form
   * <code>object:interface</code> or
   * <code>object:vector-interface[pos]</code>. Return the interface
   * name (without the <code>[pos]</code>).
   */
  static string getInterfaceFromNoun(string noun);

  /**
   * Get interface index. Decompose a string of the form
   * <code>object:interface</code> or
   * <code>object:vector-interface[pos]</code>. Return the
   * <code>pos</code> part or empty string if not present.
   */
  static string getPosArgFromNoun(string noun);

  /**
   * Return a list of the interfaces which do not have their default
   * values for the given objects.
   */
  template <typename Cont>
  static vector< pair<IBPtr, const InterfaceBase *> >
  getNonDefaultInterfaces(const Cont &);

  //@}

  /** @name Manipulate objects in the repository. */
  //@{
  /**
   * Call the InterfacedBase::update() function of all objects.
   */
  static void update();

  /**
   * Clear the InterfacedBase::touched() flag in all objects in the
   * given container.
   */
  template<typename Cont>
  static void clearAll(const Cont & c) 
  {  
    for_each(c, mem_fun(&InterfacedBase::clear));
  }

  /**
   * Set the status of all objects in the given container to
   * InterfacedBase::uninitialized.
   */
  template<typename Cont>
  static void resetAll(const Cont & c) 
  {  
    for_each(c, mem_fun(&InterfacedBase::reset));
  }

  /**
   * Setup an object. Execute the InterfacedBase::readSetup() method
   * of \a ip with the stream \a is as argument.
   */
  static void readSetup(tIBPtr ip, istream & is);

  /**
   * Lock the given object. Locked objects cannot be
   * changed through an interface.
   */
  static void lock(tIBPtr ip) { ip->lock(); }

  /**
   * Unlock the given object. Locked objects cannot be changed through
   * an interface.
   */
  static void unlock(tIBPtr ip) { ip->unlock(); }
  //@}

  /** @name Access the documentation of objects. */
  //@{
  /**
   * Return the class documentation of a given object
   */
  static const ClassDocumentationBase * getDocumentation(tcIBPtr ip);

  /**
   * Get the description for the model implemented in the class of the
   * given object.
   */
  static string getModelDescription(tcIBPtr ip);

  /**
   * Get the references for the model implemented in the class of the
   * given object.
   */
  static string getModelReferences(tcIBPtr ip);
  //@}

  /** @name Manipulate the output streams of the repository. */
  //@{
  /**
   * Set the standard output stream
   */
  static void cout(ostream & os) { coutp() = &os; }

  /**
   * Get the standard output stream
   */
  static ostream & cout() { return *coutp(); }

  /**
   * Set the standard error stream
   */
  static void cerr(ostream & os) { cerrp() = &os; }

  /**
   * Get the standard error stream
   */
  static ostream & cerr() { return *cerrp(); }

  /**
   * Set the standard log stream
   */
  static void clog(ostream & os) { clogp() = &os; }

  /**
   * Get the standard log stream
   */
  static ostream & clog() { return *clogp(); }
  //@}

protected:

  /** @name Access standard InterfacedBase functions. */
  //@{
  /**
   * Return a clone of the given object. Calls the
   * InterfacedBase::clone() function of \a t and casts the resulting
   * pointer to the correct type.
   */
  template <typename T>
  static typename Ptr<T>::pointer clone(const T & t);

  /**
   * Return a clone of the given object. Calls the
   * InterfacedBase::fullclone() function of \a t and casts the
   * resulting pointer to the correct type.
   */
  template <typename T>
  static typename Ptr<T>::pointer fullclone(const T & t);

  /**
   * Rebind references. For all objects directly referenced by \a obj,
   * replace them with the translation found in \a trans. If \a obj has a
   * Reference or a member of a RefVector interface which is null, and
   * the corresponding interface has the RefInterfaceBase::defaultIfNull() flag set,
   * translate the null pointer to the first acceptable object in
   * defaults.
   */
  static void rebind(InterfacedBase & obj, const TranslationMap & trans,
		     const IVector & defaults);
  //@}
  

  /**
   * Add interfaces to the given map for the class with the given
   * class description. Recursively do the same with the base classes.
   */
  static void addInterfaces(const ClassDescriptionBase &,
			    InterfaceMap &, bool all = true);

  /** @name Functions containing the static instances of objects used
      by the repository. */
  //@{
  /**
   * All InterfacedBase objects mapped to their name.
   */
  static ObjectMap & objects();

  /**
   * All InterfacedBase objects.
   */
  static ObjectSet & allObjects();

  /**
   * Sets of InterfaceBase objects mapped to the class description of
   * the class for which they are defined.
   */
  static TypeInterfaceMap & interfaces();

  /**
   * Sets of ClassDocumentationBase objects mapped to the class
   * description of the class for which they are defined.
   */
  static TypeDocumentationMap & documentations();

  /**
   * All defined directories.
   */
  static DirectorySet & directories();

  /**
   * The current directory stack.
   */
  static StringVector & directoryStack();

  /**
   * Flag to say if we are in the middle of an update procedure.
   */
  static bool & updating();

  /**
   * The current current standard output stream.
   */
  static ostream *& coutp();

  /**
   * The current current standard error stream.
   */
  static ostream *& cerrp();
  /**
   * The current current standard log stream.
   */
  static ostream *& clogp();
  //@}

};


}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "BaseRepository.tcc"
#endif

#endif /* ThePEG_BaseRepository_H */
