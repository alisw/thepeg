// -*- C++ -*-
#ifndef Ariadne5_ColourIndex_H
#define Ariadne5_ColourIndex_H

#include "Ariadne/Config/Ariadne5.h"

//
// This is the declaration of the ColourIndex class.
//

namespace Ariadne5 {

/**
 * The olourIndex class keeps track of the colour structure in a
 * DipoleState such that when colour reconnections are included only
 * those reconnections which are allowed can be performed.
 */
class ColourIndex {

public:

  /**
   * The default constructor.
   */
  ColourIndex(unsigned isys = 0): index(0), systemIndex(isys), parity(false) {}

public:

  /**
   * Check for equality.
   */
  bool operator==(const ColourIndex & c) const {
    return index == c.index && systemIndex == c.systemIndex;
  }

  /**
   * Check for inequality.
   */
  bool operator!=(const ColourIndex & c) const {
    return index != c.index || systemIndex != c.systemIndex;
  }

  /**
   * Check for ordering.
   */
  bool operator<(const ColourIndex & c) const {
    return index < c.index || ( index == c.index && systemIndex < c.systemIndex );
  }

  /**
   * Return the system index.
   */
  unsigned system() const {
    return systemIndex;
  }

  /**
   * Set the system index.
   */
  void system(unsigned i) {
    systemIndex = i;
  }

  /**
   * Randomize this index making sure it is not the same as any of the
   * indices given as argument.
   */
  void generate(const ColourIndex & c1 = ColourIndex(),
		const ColourIndex & c2 = ColourIndex(),
		const ColourIndex & c3 = ColourIndex());


  /**
   * Generate a new index making sure it is not the same as this index
   * or any of the indices given as argument.
   */
  ColourIndex generateNew(const ColourIndex & c1 = ColourIndex(),
			  const ColourIndex & c2 = ColourIndex()) const {
    ColourIndex ret(system());
    ret.generate(*this, c1, c2);
    return ret;
  }

  /**
   * Print function for PersistentOStream.
   */
  template <typename OS>
  void print(OS & os) const {
    os << index << systemIndex << parity;
  }

  /**
   * Read function for PersistentIStream.
   */
  template <typename IS>
  void read(IS & is) {
    is >> index >> systemIndex >> parity;
  }



private:

  /**
   * The actual colour index.
   */
  unsigned index;

  /**
   * An additional index to keep different systems separate.
   */
  unsigned systemIndex;

  /**
   * A parity for allowing junction swings.
   */
  bool parity;

};

}

// template <typename OS>
// OS & operator<<(OS & os, const Ariadne5::ColourIndex & c) {
//   c.print(os);
//   return os;
// }

// template <typename IS>
// IS & operator>>(IS & is, Ariadne5::ColourIndex & c) {
//   c.read(is);
//   return is;
// }

inline std::ostream & operator<<(std::ostream & os, const Ariadne5::ColourIndex & c) {
  c.print(os);
  return os;
}

inline ThePEG::PersistentOStream & operator<<(ThePEG::PersistentOStream & os,
					      const Ariadne5::ColourIndex & c) {
  c.print(os);
  return os;
}

template <typename IS>
IS & operator>>(IS & is, Ariadne5::ColourIndex & c) {
  c.read(is);
  return is;
}

#endif /* Ariadne5_ColourIndex_H */
