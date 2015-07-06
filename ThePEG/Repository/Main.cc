// -*- C++ -*-
//
// Main.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Main class.
//

#include "Main.h"

using namespace ThePEG;

AbstractNoPIOClassDescription<Main> Main::initMain;
// Definition of the static class description member.

EGPtr Main::theEventGenerator;

long Main::theN = 0;

vector<string> Main::theArguments;

