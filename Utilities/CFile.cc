// -*- C++ -*-
//
// CFileLineReader.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CFileLineReader class.
//

#include "CFile.h"
#include "config.h"
#include <cstdlib>
#include <cstdio>
#include "Throw.h"
#include <cstring>

#ifdef HAVE_LIBZ
#include <zlib.h>
#else
typedef void * gzFile;
#endif
#ifdef HAVE_LIBBZ2_NEVER
#include <bzlib.h>
#else
typedef void BZFILE;
#endif

using namespace ThePEG;

void CFile::open(string filename, string mode) {
  close();
  if ( filename[filename.length()-1] == '|' &&
       mode.find("r") != string::npos ) {
    filename = filename.substr(0, filename.length() - 1);
    file = popen(filename.c_str(), mode.c_str());
    fileType = pipe;
  }
  else if ( filename[0] == '|' &&
	    mode.find("w") != string::npos ) {
    filename = filename.substr(1);
    file = popen(filename.c_str(), mode.c_str());
    fileType = pipe;
  }
  else if ( filename.substr(filename.length()-3,3) == ".gz" ) {
#ifdef HAVE_LIBZ
    file = gzopen(filename.c_str(), mode.c_str());
    fileType = gzip;
#else
#ifdef ThePEG_GZREAD_FILE
#ifdef ThePEG_GZWRITE_FILE
    if ( mode.find("r") != string::npos )
      filename = ThePEG_GZREAD_FILE " " + filename + " 2>/dev/null";
    else
      filename = ThePEG_GZWRITE_FILE " " + filename + " 2>/dev/null";
    file = popen(filename.c_str(), mode.c_str());
    fileType = pipe;
#else
    file = fopen(filename.c_str(), mode.c_str());
    fileType = plain;
#endif
#endif
#endif
  }
  else if ( filename.substr(filename.length()-4,4) == ".bz2" ) {
#ifdef HAVE_LIBBZ2_NEVER
    file = BZ2_bzopen(filename.c_str(), mode.c_str());
    fileType = bzip2;
#else
#ifdef ThePEG_BZ2READ_FILE
#ifdef ThePEG_BZ2WRITE_FILE
    if ( mode.find("r") != string::npos )
      filename = ThePEG_BZ2READ_FILE " " + filename + " 2>/dev/null";
    else
      filename = ThePEG_BZ2WRITE_FILE " " + filename + " 2>/dev/null";
    file = popen(filename.c_str(), mode.c_str());
    fileType = pipe;
#else
    file = fopen(filename.c_str(), mode.c_str());
    fileType = plain;
#endif
#endif
#endif
  }
  else {
    file = fopen(filename.c_str(), mode.c_str());
    fileType = plain;
  }
  if ( !file ) {
    Throw<FileError>() 
      << std::strerror(errno) << ": " << filename 
      << Exception::runerror;
  }
}

void CFile::close() {
  if ( !file ) {
    fileType = undefined;
    return;
  }
  switch ( fileType ) {
  case plain:
    fclose((FILE*)file);
    break;
  case pipe:
    pclose((FILE*)file);
    break;
#ifdef HAVE_LIBZ
  case gzip:
    gzclose((gzFile)file);
    break;
#endif
#ifdef HAVE_LIBBZ2_NEVER
  case bzip2:
    BZ2_bzclose(file);
    break;
#endif
  default:
    break;
  }
  file = 0;
  fileType = undefined;
}

char * CFile::gets(char * s, int size) {
  switch ( fileType ) {
  case plain:
  case pipe: return fgets(s, size, (FILE*)file);
#ifdef HAVE_LIBZ
  case gzip: return gzgets((gzFile)file, s, size);
#endif
#ifdef HAVE_LIBBZ2_NEVER
  case bzip2: // don't know what to do here
#endif
  default:
    return 0;
  }
}

int CFile::puts(const char * s) {
  switch ( fileType ) {
  case plain:
  case pipe: return fputs(s, (FILE*)file);
#ifdef HAVE_LIBZ
  case gzip: return gzputs((gzFile)file, s);
#endif
#ifdef HAVE_LIBBZ2_NEVER
  case bzip2: // don't know what to do here
#endif
  default:
    return 0;
  }
}

int CFile::getc() {
  switch ( fileType ) {
  case plain:
  case pipe: return fgetc((FILE*)file);
#ifdef HAVE_LIBZ
  case gzip: return gzgetc((gzFile)file);
#endif
#ifdef HAVE_LIBBZ2_NEVER
    case bzip2:// don't know what to do here
#endif
  default:
    return 0;
  }
}

int CFile::putc(int c) {
  switch ( fileType ) {
  case plain:
  case pipe: return fputc(c, (FILE*)file);
#ifdef HAVE_LIBZ
  case gzip: return gzputc((gzFile)file, c);
#endif
#ifdef HAVE_LIBBZ2_NEVER
  case bzip2: // don't know what to do here
#endif
  default:
    return 0;
  }
}

int CFile::ungetc(int c) {
  switch ( fileType ) {
  case plain:
  case pipe: return std::ungetc(c, (FILE*)file);
#ifdef HAVE_LIBZ
  case gzip: return gzungetc(c, (gzFile)file);
#endif
#ifdef HAVE_LIBBZ2_NEVER
  case bzip2: // don't know what to do here
#endif
  default:
    return 0;
  }
}

size_t CFile::read(void *ptr, size_t size, size_t nmemb) {
  switch ( fileType ) {
  case plain:
  case pipe: return fread(ptr, size, nmemb, (FILE*)file);
#ifdef HAVE_LIBZ
  case gzip: return gzread((gzFile)file, ptr, size);
#endif
#ifdef HAVE_LIBBZ2_NEVER
  case bzip2: return BZ2_bzread(file, ptr, size);
#endif
  default:
    return 0;
  }
}

size_t CFile::write(const void *ptr, size_t size, size_t nmemb) {
  switch ( fileType ) {
  case plain:
  case pipe: return fwrite(ptr, size, nmemb, (FILE*)file);
#ifdef HAVE_LIBZ
  case gzip:  return gzwrite((gzFile)file, ptr, size);
#endif
#ifdef HAVE_LIBBZ2_NEVER
  case bzip2: return BZ2_bzwrite(file, ptr, size);
#endif
  default:
    return 0;
  }
}



