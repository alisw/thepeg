// -*- C++ -*-
#ifndef FragConfig_H
#define FragConfig_H

#include "Pythia7/Config/Pythia7.h"
#include "ThePEG/Vectors/LorentzVector.h"
#include "ThePEG/Vectors/Transverse.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Config/Pointers.h"

namespace Pythia7{ 

class LundFragmentationHandler;
class EndPoint;

class StringRegion;
/** Bare pointer to StringRegion. */
typedef StringRegion* StringRegionPtr;
/** Bare pointer to const StringRegion. */
typedef const StringRegion* cStringRegionPtr; 

class String;
/** Bare pointer to String. */
typedef String* StringPtr; 
/** Bare pointer to const String. */
typedef const String* cStringPtr; 

/** Iterator in a list of Particle pointers. */
typedef ParticleList::iterator ParticleListIt; 
/** Const iterator in a list of Particle pointers. */
typedef ParticleList::const_iterator cParticleListIt; 

/** Iterator in a list of transient Particle pointers. */
typedef tPList::iterator tPListIt; 
/** Const iterator in a list of transient Particle pointers. */
typedef tPList::const_iterator tPListCIt; 

ThePEG_DECLARE_CLASS_POINTERS(LundPtGenerator, LPtGenPtr);
ThePEG_DECLARE_CLASS_POINTERS(LundZGenerator, LZGenPtr);
ThePEG_DECLARE_CLASS_POINTERS(LundFragHandler, LFragHdlrPtr);
ThePEG_DECLARE_CLASS_POINTERS(LundFlavourGenerator, LFlGenPtr);
ThePEG_DECLARE_CLASS_POINTERS(LundFlavourGenerator, FlavourGeneratorPtr);


/** A vector of four momenta. */
typedef vector<LorentzMomentum> MomentumVector;

using namespace ThePEG::Math;


}//Pythia7EndNamespace

namespace ThePEG {

ThePEG_DECLARE_CLASS_POINTERS(PtGenerator, PtGeneratorPtr);
ThePEG_DECLARE_CLASS_POINTERS(ZGenerator, ZGeneratorPtr);

}

#endif // FragConfig_H 
