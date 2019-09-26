// -*- C++ -*-
#ifndef ThePEG_BlobMEBase_H
#define ThePEG_BlobMEBase_H
//
// This is the declaration of the BlobMEBase class.
//

#include "ThePEG/MatrixElement/MEBase.h"

namespace ThePEG {

/**
 * Here is the documentation of the BlobMEBase class.
 *
 * @see \ref BlobMEBaseInterfaces "The interfaces"
 * defined for BlobMEBase.
 */
class BlobMEBase: public MEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BlobMEBase();

  /**
   * The destructor.
   */
  virtual ~BlobMEBase();
  //@}

public:

  /**
   * Helper struct to represent colour connections.
   */
  struct ColourConnection {

    /**
     * The members of the colour connection
     */
    vector<int> members;

    /**
     * Add a leg's colour to the connection
     */
    void addColour(int leg) {
      members.push_back(leg+1);
    }

    /**
     * Add a leg's anti-colour to the connection
     */
    void addAntiColour(int leg) {
      members.push_back(-leg-1);
    }

    /**
     * Write out the connection to the colour lines string
     */
    string write(size_t& sourceCount, bool sink) const;

  };

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;
  //@}

  /**
   * Return the possible processes this matrix element will be able to handle,
   * as a map incoming to outgoing; it is assumed that the number of outgoing
   * partons does not vary.
   */
  virtual multimap<tcPDPair,tcPDVector> processes() const = 0;

  /**
   * Return the colour connections for the process.
   */
  virtual list<ColourConnection> colourConnections() const = 0;

public:

  /**
   * Describe an abstract base class without persistent data.
   */
  static AbstractNoPIOClassDescription<BlobMEBase> initBlobMEBase;

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BlobMEBase & operator=(const BlobMEBase &) = delete;

  /**
   * The colour lines object used as a proxy to connect colours in
   * BlobDiagram::construct
   */
  mutable ColourLines theColourLines;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * BlobMEBase.
 */
template <>
struct BaseClassTrait<BlobMEBase,1>: public ClassTraitsType {
  /** Typedef of the base class of BlobMEBase. */
  typedef MEBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * BlobMEBase class.
 */
template <>
struct ClassTraits<BlobMEBase>: public ClassTraitsBase<BlobMEBase> {
  /** Return the class name. */
  static string className() { return "ThePEG::BlobMEBase"; }
};

/** @endcond */

}


#endif /* ThePEG_BlobMEBase_H */
