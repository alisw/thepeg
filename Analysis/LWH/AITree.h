// -*- C++ -*-
//
// AITree.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AITree_H
#define LWH_AITree_H



/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IManagedObject;

class ITree {

public:

  virtual ~ITree() {}

  virtual std::string storeName() const = 0;
  virtual IManagedObject * find(const std::string & name) = 0;
  virtual std::string pwd() const = 0;
  virtual bool commit() = 0;
  virtual bool close() = 0;
  virtual bool mkdir(const std::string &) = 0;
  virtual bool mkdirs(const std::string &) = 0;
  virtual bool cd(const std::string &) = 0;
  virtual bool rmdir(const std::string & str) = 0;
  virtual bool rm(const std::string & str) = 0;
  virtual std::string findPath(const IManagedObject & o) const = 0;
  virtual bool mv(const std::string & oldo, const std::string & newo) = 0;
  virtual void setOverwrite(bool o = true) = 0;
  virtual bool cp(const std::string &, const std::string &, bool = false) = 0;

};

}

/** @endcond */





#endif /* LWH_AITree_H */
