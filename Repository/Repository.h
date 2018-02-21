// -*- C++ -*-
//
// Repository.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Repository_H
#define ThePEG_Repository_H
// This is the declaration of the Repository class.

#include "ThePEG/Config/ThePEG.h"
#include "BaseRepository.h"
#include "EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/MatcherBase.h"

namespace ThePEG {

/**
 * Repository inherits from the BaseRepository class. While
 * BaseRepository is fairly general and could in principle be used for
 * any program where sets of InterfacedBase objects are managed, the
 * Repository is ThePEG specific in that it deals with ParticleData,
 * ParticleMatchers and EventGenerators.
 *
 * One main function is to write an EventGenerator to disk using
 * saveRun(). Here all objects needed for the run, including the
 * EventGenerator is cloned and isolated from the other objects in the
 * Repository (and are hence not handled by the Repository anymore)
 * before they are all persistently written out to disk.
 *
 * The Register() function simply pass the objects to the corresonding
 * method in BaseRepository, but if the object is a ParticleData or a
 * ParticleMatcher, they are stored separately.
 *
 * @see BaseRepository
 * @see InterfacedBase
 * @see ParticleData
 * @see ParticleMatcher
 * @see EventGenerator
 * 
 */
class Repository: public BaseRepository {

public:

  /** A map of EventGenerator objects indexed by their run name. */
  typedef map<string,EGPtr> GeneratorMap;

public:

  /** @name Standsrd constructors and destructors */
  //@{
  /**
   * The default constructor is the only one that should be used.
   */
  Repository();

  /**
   * The destructor will do some clean-up when the last Repository is
   * deleted.
   */
  ~Repository();

public:

  /** @name Functions for register objects in the Repository. */
  //@{
  /**
   * Register an object with BaseRepository::Register() and add it to
   * the list of particles or matchers if of any of those
   * types.
   */
  static void Register(IBPtr);

  /**
   * Register an object with BaseRepository::Register() and add it to
   * the list of particles or matchers if of any of those
   * types.
   */
  static void Register(IBPtr, string newName);
  //@}

  /** @name Access ParticleData and MatcherBase objects in the
      Repository. */
  //@{
  /**
   * Add a particle to the list of default ones. If one of the same
   * type alredy existed, it is removed from the list (but not from
   * the repository).
   */
  static void defaultParticle(tPDPtr);

  /**
   * Get a pointer to the default particle of the given type or
   * generic name.
   */
  static PDPtr defaultParticle(PID id);

  /**
   * Get a pointer to a particle based on the given path or name. The
   * argument is first treated as a path to an object. If no such
   * particle object is found, the argument is treated as a generic
   * particle PDGName and is searched for among the default set of
   * particles.
   */
  static tPDPtr findParticle(string name);

  /**
   * Return the set of all particles in the repository.
   */
  static const ParticleDataSet & allParticles() { return particles(); }

  /**
   * Return the set of all matchers in the repository.
   */
  static const MatcherSet & allMatchers() { return matchers(); }

  /**
   * Find a matcher with a given generic name
   */
  static tPMPtr findMatcher(string name);

  /**
   * Special function for copying particles. Also corresponding
   * anti-particle is copied to the same directory. In addition, their
   * decay modes are copied.
   */
  static string copyParticle(tPDPtr, string);
  //@}

  /** @name Functions to isolate Eventgenerator objects. */
  //@{
  /**
   * Isolate an event generator, \a eg, and save it to disk in a file
   * named \a name (with <code>.run</code> appended.
   */
  static EGPtr makeRun(tEGPtr eg, string name);

  /**
   * Isolate an event generatorn, named \a EGname, set its run \a name
   * and save it to a file named \a filename.
   */
  static void saveRun(string EGname, string name, string filename);
  //@}

  /** @name I/O functions for the Repository. */
  //@{
  /**
   * Load a whole repository from the given file. All objects
   * previously in the Repository are discarded. Any errors will be
   * reported in the returned string.
   */
  static string load(string filename);

  /**
   * Save the repository to the given file.
   */
  static void save(string filename);

  /**
   * Save the repository to the default file.
   */
  static void save() { save(currentFileName()); }

  /**
   * Write some statistics about the repository to the standard output.
   */
  static void stats(ostream &);
  //@}

  /** @name Command-line interface functions. */
  //@{
  /**
   * Print out a help message. Extended text for a specific command if given.
   */
  static void help(string command, ostream & os);

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
   * Read commands from a stream and send them one by one to exec().
   *
   * @param is the stream from which to read commands.
   * @param os the stream where output is written.
   * @param prompt before reading a command from \a is, this string is
   * written to \a os.
   */
  static void read(istream & is, ostream & os, string prompt = "");

  /**
   * Read commands from a file and send them one by one to exec().
   *
   * Passes the call through to read(istream, ostream), but also sets
   * currentReadDirStack() correctly.
   *
   * Returns possible messages.
   *
   * @param filename the file from which to read commands.
   * @param os the stream where output is written.
   */
  static string read(string filename, ostream & os);

  /**
   * Interpret the command in \a cmd and return possible
   * messages. This is the main function for the command-line
   * interface. The syntax is described elsewhere. The ostream
   * argument is currently unused.
   */
  static string exec(string cmd, ostream &);

  /**
   * Insert the given EventGenerator and its dependent Interfaced
   * objects into the repository and read commands to modify its
   * interfaces. Any line accepted by the command-line interface will
   * be executed, but the main purpose of this function is to modify
   * an already saved and initialized EventGenerator before running
   * without re-initializing. If an interface which does not have the
   * dependencySafe() flag set, a warning will be emitted.
   */
  static string modifyEventGenerator(EventGenerator & eg, string filename, 
				     ostream & os, bool initOnly = false);  

  /**
   * Reset the given EventGenerator; equivalent to
   * modifyEventGenerator without reading an input file.
   */
  static void resetEventGenerator(EventGenerator & eg);  
  //@}

  /**
   * Return the version number of ThePEG.
   */
  static string version();

  /**
   * Return a string with a ThePEG banner.
   */
  static string banner();

private:

  /**
   * Used by Register.
   */
  static void registerParticle(tPDPtr);

  /**
   * Used by Register.
   */
  static void registerMatcher(tPMPtr);

  /** 
   * Used by read()
   */
  static void execAndCheckReply(string, ostream &);

  /**
   *  Check that the PDG name is not a duplicate
   */
  static void checkDuplicatePDGName(PDPtr);

protected:

  /** @name Functions containing the static instances of objects used
      by the repository. */
  //@{
  /**
   * The set of default particles.
   */
  static ParticleMap & defaultParticles();

  /**
   * The set of all particles.
   */
  static ParticleDataSet & particles();

  /**
   * The set of all matchers.
   */
  static MatcherSet & matchers();

  /**
   * All isolated generators mapped to their run name.
   */
  static GeneratorMap & generators();

  /**
   * The default file name used by save().
   */
  static string & currentFileName();

public:

  /**
   * If non-zero the setup program will exit with this error code as
   * soon as an error is encountered.
   */
  static int & exitOnError();

  /**
   * Call this function to clean up the repository at the end of your
   * program if you are using the static functions directly without
   * going through a Repository object. There, the destructor would do
   * the job.
   */
  static void cleanup();
  //@}

private:

  /**
   * It makes no sense to copy a Repository, so this constructor is
   * not implemented
   */
  Repository(const Repository &);

  /**
   * It makes no sense to copy a Repository, so this assignment is
   * not implemented
   */
  Repository & operator=(const Repository &);

  /**
   * Count the number of repositorys instantiated.
   */
  static int ninstances;
  

};

}

#endif /* ThePEG_Repository_H */
