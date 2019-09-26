// -*- C++ -*-
//
// VSSVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VSSVertex_H
#define ThePEG_VSSVertex_H
//
// This is the declaration of the VSSVertex class.
//
#include "GeneralVSSVertex.h"
#include "VSSVertex.fh"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *
 *  The VSSVertex class is the implementation of the vector-scalar-scalar
 *  vertex. It inherits from the AbstractVSSVertex class for storage of the particles 
 *  and implements the helicity calculations.
 *
 *  All such vertices should inherit from this class and implement the virtual
 *  setCoupling member
 *
 *  The form of the vertex is
 * \f[-ic\left(p_2-p_3\right)\cdot\epsilon_1\phi_2\phi_3\f]
 *
 *  @see GeneralVSSVertex
 */
class VSSVertex: public GeneralVSSVertex {
    
public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

private:
  
  /**
   * Private and non-existent assignment operator.
   */
  VSSVertex & operator=(const VSSVertex &) = delete;
  
};

}

}
#endif /* ThePEG_VSSVertex_H */
