// -*- C++ -*-
//
// LesHouchesFileReader.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_LesHouchesFileReader_H
#define THEPEG_LesHouchesFileReader_H
// This is the declaration of the LesHouchesFileReader class.

#include "LesHouchesReader.h"
#include "LesHouchesFileReader.fh"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/Utilities/CFileLineReader.h"
#include <string>
#include <stdio.h>

namespace ThePEG {


/**
 * LesHouchesFileReader is an base class to be used for objects which
 * reads event files from matrix element generators. It inherits from
 * LesHouchesReader and extends it by defining a file handle to be
 * read from, which is opened and closed by the open() and close()
 * functions. Note that the file handle is a standard C filehandle and
 * not a C++ stream. This is because there is no standard way in C++
 * to connect a pipe to a stream for reading eg. gzipped files. This
 * class is able to read plain event files conforming to the Les
 * Houches Event File accord.
 *
 * @see \ref LesHouchesFileReaderInterfaces "Th1e interfaces"
 * defined for LesHouchesFileReader.
 * @see Event
 * @see LesHouchesReader
 */
class LesHouchesFileReader: public LesHouchesReader {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  LesHouchesFileReader() : neve(0), ieve(0), theQNumbers(false) {}

  /**
   * Copy-constructor. Note that a file which is opened in the object
   * copied from will have to be reopened in this.
   */
  LesHouchesFileReader(const LesHouchesFileReader &);

  /**
   * Destructor.
   */
  virtual ~LesHouchesFileReader();
  //@}

public:

  /** @name Virtual functions specified by the LesHouchesReader base class. */
  //@{
  /**
   * Initialize. This function is called by the LesHouchesEventHandler
   * to which this object is assigned.
   */
  virtual void initialize(LesHouchesEventHandler & eh);

  /**
   * Open a file with events. Derived classes should overwrite it and
   * first calling it before reading in the run information into the
   * corresponding protected variables.
   */
  virtual void open();

  /**
   * Close the file from which events have been read.
   */
  virtual void close();
 

  /**
   * Read the next event from the file or stream into the
   * corresponding protected variables. Return false if there is no
   * more events or if this was not a LHF event file.
   */
  virtual bool doReadEvent();
  //@}

  /**
   * Return the name of the file from where to read events.
   */
  string filename() const { return theFileName; }

  /** 
   * Return the optional weights information string ("Names")
   */

  virtual vector<string> optWeightsNamesFunc();
  

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();


  /** 
   * Erases all occurences of a substring from a string 
   */
  
  void erase_substr(std::string& subject, const std::string& search);


protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

  /** @name Standard (and non-standard) Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Return true if this object needs to be initialized before all
   * other objects because it needs to extract PDFs from the event file.
   */
  virtual bool preInitialize() const;
  //@

protected:

  /**
   * The wrapper around the C FILE stream from which to read
   */
  CFileLineReader cfile;

protected:

  /**
   * The number of events in this file.
   */
  long neve;

  /**
   * The current event number.
   */
  long ieve;

  /**
   * If the file is a standard Les Houches formatted file (LHF) this
   * is its version number. If empty, this is not a Les Houches
   * formatted file
   */
  string LHFVersion;

  /**
   * If LHF. All lines (since the last open() or readEvent()) outside
   * the header, init and event tags.
   */
  string outsideBlock;

  /**
   * If LHF. All lines from the header block.
   */
  string headerBlock;

  /**
   * If LHF. Additional comments found in the init block.
   */
  string initComments;

  /**
   * If LHF. Map of attributes (name-value pairs) found in the init
   * tag.
   */
  map<string,string> initAttributes;

  /**
   * If LHF. Additional comments found with the last read event.
   */
  string eventComments;

  /**
   * If LHF. Map of attributes (name-value pairs) found in the last
   * event tag.
   */
  map<string,string> eventAttributes;

private:

  /**
   * The name of the file from where to read events.
   */
  string theFileName;

  /**
   *  Whether or not to search for QNUMBERS stuff
   */
  bool theQNumbers;

  /**
   *  Decayer for any decay modes read from the file
   */
  DecayerPtr theDecayer;
  
  /**
   * Further information on the weights
   */
  map<string,string> scalemap;

  /**
   * Temporary holder for optional weights
   */
  
  map<string,double> optionalWeightsTemp;


private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static ClassDescription<LesHouchesFileReader> initLesHouchesFileReader;

  /**
   * Private and non-existent assignment operator.
   */
  LesHouchesFileReader & operator=(const LesHouchesFileReader &);

public:

  /** @cond EXCEPTIONCLASSES */
  /** Exception class used by LesHouchesFileReader if reading the file
   *  fails. */
  class LesHouchesFileError: public Exception {};
  /** @endcond */

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of LesHouchesFileReader.
 */
template <>
struct BaseClassTrait<LesHouchesFileReader,1>: public ClassTraitsType {
  /** Typedef of the base class of LesHouchesFileReader. */
  typedef LesHouchesReader NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * LesHouchesFileReader class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<LesHouchesFileReader>
  : public ClassTraitsBase<LesHouchesFileReader> {
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::LesHouchesFileReader"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the LesHouchesFileReader class and every other class it uses
   * (except the base class).
   */
  static string library() { return "LesHouches.so"; }

};

/** @endcond */

}

#endif /* THEPEG_LesHouchesFileReader_H */
