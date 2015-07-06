// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Emission class.
//

#include "Emission.h"
#include "EmitterBase.h"
#include "Parton.h"
#include "DipoleState.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/UnitIO.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "Ariadne/Config/UnitFO.h"

using namespace Ariadne5;

void Emission::fillReferences(CloneSet & cset) const {}

void Emission::rebind(const TranslationMap & trans) {
  dipole = trans.translate(dipole);
  cdipole = trans.translate(cdipole);
  vector<tParPtr> rold;
  rold.swap(radiators);
  trans.translate(inserter(radiators), rold.begin(), rold.end());
  vector<tParPtr> parold;
  parold.swap(partons);
  trans.translate(inserter(partons), parold.begin(), parold.end());
  vector<tParPtr> aold;
  aold.swap(affected);
  trans.translate(inserter(affected), aold.begin(), aold.end());
  colourParent = trans.translate(colourParent);
  antiColourParent = trans.translate(antiColourParent);
  mainParent = trans.translate(mainParent);
}

bool Emission::perform() const {
  state = performed;
  if ( !model->perform(*this) ) return false;
  if ( model->finalVeto(*this) ) return false;
  return true;
}

void Emission::revert() const {
  
  model->revert(*this);
  state = reverted;
}

void Emission::persistentOutput(PersistentOStream & os) const {
  os << oenum(state) << emno << geno
     << ounit(pold.first, GeV) << ounit(pold.second, GeV)
     << ounit(rho, GeV) << y << ounit(genmom.first, GeV)
     << ounit(genmom.second, GeV) << ounit(genmom.third, GeV)
     << ymax << orderAlphaS << orderAlphaEW
     << dipole << cdipole << radiators << partons << affected << prob
     << weightPDF << reversible << failsafe << colourParent << antiColourParent
     << mainParent;
}

void Emission::persistentInput(PersistentIStream & is, int) {
  is >> ienum(state) >> emno >> geno
     >> iunit(pold.first, GeV) >> iunit(pold.second, GeV)
     >> iunit(rho, GeV) >> y  >> iunit(genmom.first, GeV)
     >> iunit(genmom.second, GeV) >> iunit(genmom.third, GeV)
     >> ymax >> orderAlphaS >> orderAlphaEW
     >> dipole >> cdipole >> radiators >> partons >> affected >> prob
     >> weightPDF >> reversible >> failsafe >> colourParent >> antiColourParent
     >> mainParent;
}

DescribeAbstractClass<Emission,CloneBase>
describeAriadne5Emission("Ariadne5::Emission", "libAriadne5.so");

void Emission::Init() {}

void Emission::debugme() const {
  LorentzMomentum diff = genmom.first + genmom.second + genmom.third
    - pold.first - pold.second;
  cerr << endl
       << "Emission type: " << StringUtils::typeName(typeid(*this)) << endl
       << " Emitter type: " << StringUtils::typeName(typeid(*model)) << endl
       << "Generation no: " << geno << endl
       << "  Emission no: " << emno << endl
       << "        State: ";
  switch ( state ) {
  case generated:
    cerr << "generated";
    break;
  case performed:
    cerr << "performed";
    break;
  case reverted:
    cerr << "reverted";
  }
  if ( failsafe ) cerr << " (failsafe)";
  cerr << endl
       << "        scale: " << rho/GeV << endl
       << "         pold: " << founit(pold.first, GeV, 10, 4) << endl
       << "               " << founit(pold.second, GeV, 10, 4) << endl
       << "       genmom: " << founit(genmom.first, GeV, 10, 4) << endl
       << "               " << founit(genmom.second, GeV, 10, 4) << endl
       << "               " << founit(genmom.third, GeV, 10, 4) << endl
       << "         diff: " << founit(diff, GeV, 10, 1, ios::scientific) << endl
       << "       dipole: " << cdipole->state()->index(cdipole) << endl
       << "      partons:";
  for ( unsigned i = 0; i < radiators.size(); ++i )
    cerr << (radiators[i] == mainParent? " *": "  ")
	 << cdipole->state()->index(radiators[i]);
  cerr << " ->";
  for ( unsigned i = 0; i < partons.size(); ++i )
    cerr << " " << cdipole->state()->index(partons[i]);
  cerr << endl;

}
