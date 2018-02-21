// -*- C++ -*-
//
// TreeFactory.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_TreeFactory_H
#define LWH_TreeFactory_H
//
// This is the declaration of the TreeFactory class.
//

#include "AITreeFactory.h"
#include <string>
#include <stdexcept>
#include "Tree.h"

namespace LWH {

using namespace AIDA;

/**
 * The creator of trees.
 */
class TreeFactory: public ITreeFactory {

public:

  /// Destructor.
  virtual ~TreeFactory() {
    clear();
  }

  /**
   * Creates a new tree that is not associated with a store.
   */
  ITree * create() {
    Tree * tree = new Tree;
    trees.insert(tree);
    return tree;
  }

  /**
   * Creates a new Tree and associates it with a store.
   * The store is assumed to be write-only.
   * The store will be created.
   * @param storeName The name of the store, if empty (""), the tree is
   *                  created in memory and therefore will not be associated
   *                  with a file.
   */
  Tree * createTree(const std::string & storeName) {
    return new Tree(storeName);
  }

  /**
   * Creates a new Tree and associates it with a store.
   * The store is assumed to be write-only.
   * The store will be created.
   * @param storeName The name of the store, if empty (""), the tree is
   *                  created in memory and therefore will not be associated
   *                  with a file.
   * @param storeType must be "xml".
   * @param readOnly  must be false since we cannot read in trees.
   * @param createNew must be true indicating that the file will be created
   */
  ITree * create(const std::string & storeName,
		 const std::string & storeType = "",
		 bool readOnly = false, bool createNew = false,
		 const std::string & = "") {
    if ( storeType != "xml" && storeType != "" && storeType != "flat" )
      throw std::runtime_error("Can only store trees in xml or flat format.");
    if ( readOnly || !createNew )
      throw std::runtime_error("Cannot read in trees.");
    return new Tree(storeName, storeType != "flat");
  }

private:

  /** Delete all trees. */
  void clear() {
    for ( std::set<Tree *>::iterator it = trees.begin();
	  it != trees.end(); ++it ) delete *it;
    trees.clear();
  }

  /** The created trees. */
  std::set<Tree *> trees;

};

}

#endif /* LWH_TreeFactory_H */
