// -*- C++ -*-
#ifndef THEPEG_CFile_H
#define THEPEG_CFile_H
//
// This is the declaration of the CFile class.
//

#include "Exception.h"

namespace ThePEG {

/**
 * Here is the documentation of the CFile class.
 */
class CFile {

public:

  /**
   *  Type of the file
   */
  enum FileType {
    undefined, plain, pipe, gzip, bzip2
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CFile(): file(0), fileType(undefined) {}

  /**
   * Create a CFile given a file name and a mode.
   */
  CFile(string filename, string mode)
  : file(0), fileType(undefined) {
    open(filename, mode);
  }
    
  /**
   * The destructor.
   */
  ~CFile() {}
  //@}

  /**
   * Open the file
   */
  void open(string filename, string mode);

  /**
   *  Close the file
   */
  void close();

  /**
   *  Pointer to the file
   */
  operator void * () const {
    return fileType != undefined? file: 0;
  }

  /**
   *  Exist for file existance
   */
  bool operator!() const {
    return !(operator void * ());
  }

  /**
   *  Get characters
   */
  char * gets(char * s, int size);

  /**
   *  Set characters
   */
  int puts(const char * s);

  /**
   *  Get the current character
   */
  int getc();

  /**
   *  Set the current character
   */
  int putc(int c);

  /** Pushes the byte specified by \a c (converted to an unsigned char)
   *  back onto the stream
   */
  int ungetc(int c);

  /**
   *  Read
   */
  size_t read(void *ptr, size_t size, size_t nmemb = 1);

  /**
   *   Write
   */
  size_t write(const void *ptr, size_t size, size_t nmemb = 1);

private:

  /**
   * Pointer to the file
   */
  void * file;

  /**
   *  Type of the file
   */
  FileType fileType;

public:

  /** @cond EXCEPTIONCLASSES */
  /** Exception class used by CFile if reading a file fails. */
  class FileError: public Exception {};
  /** @endcond */

};

}


#endif /* THEPEG_CFile_H */
