// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CloneBase class.
//

#include "CloneBase.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Ariadne5;

CloneBase::~CloneBase() {}

void CloneBase::fillReferences(CloneSet &) const {}

void CloneBase::rebind(const TranslationMap & trans) {}

// Definition of the static class description member.
DescribeAbstractNoPIOClass<CloneBase,PersistentBase>
describeAriadne5CloneBase("Ariadne5::CloneBase", "libAriadne5.so");

void CloneBase::Init() {}
