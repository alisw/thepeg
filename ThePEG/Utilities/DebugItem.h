// -*- C++ -*-
#ifndef ThePEG_DebugItem_H
#define ThePEG_DebugItem_H
//
// This is the declaration of the DebugItem class.
//

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * The DebugItem class can be used to efficiently handle detailed
 * debug options. The actual objects are used anywhere in a function
 * where optional debugging should be done. At that point a static
 * object of DebugItem should be constructed giving a name to be used
 * (it should be static to ensure that the initialization is only done
 * once). After that the object is automatically cast to a bool
 * indicating whether or not debugging has been requested for this
 * item.
 */
class DebugItem {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The only relevant constructor. The string should typically be on
   * the form "namspace::class::function". If the Debug::level is
   * larger than or equal to the given \a level this DebugItem will be
   * turned on.
   */
  DebugItem(string itemname, int level = 100);
  //@}

public:

  /**
   * Switch on all DebugItem objects matching the given string. If \a
   * after is positive delay the DebugItem until that number of
   * tics. If the string is on the form "<name>=<int>" the integer
   * will be taken as the delay.
   */
  static void setDebugItem(string itemname, long after = 0);

  /**
   * Advance one tic, opssibly switching on more debug items.
   */
  static void tic();

  /**
   * Cheap way of testing if debugging should be done.
   */
  operator bool () const {
#ifndef ThePEG_NO_DEBUG
    return debug;
#else
    return false;
#endif
  }

private:

  /**
   * Set to true if debugging requested.
   */
  bool debug;

  /**
   * Counter for number of tics.
   */
  static long & ticker();

  /**
   * The DebugItem objects registered, indexed by their name.
   */
  static multimap<string,DebugItem*> & items();

  /**
   * The DebugItem objects registered, indexed by the tic at which
   * they should be turned on..
   */
  static multimap<long,DebugItem*> & itemtics();

  /**
   * The DebugItem names registered together with the tic at which it
   * should be turned on.
   */
  static map<string,long> & nametics();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DebugItem & operator=(const DebugItem &);

};

}

#endif /* ThePEG_DebugItem_H */
