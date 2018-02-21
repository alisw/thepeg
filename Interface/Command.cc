// -*- C++ -*-
//
// Command.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CommandBase class.
//

#include "InterfacedBase.h"
#include "Command.h"

using namespace ThePEG;

string CommandBase::exec(InterfacedBase & i, string action,
			 string arguments) const {
  if ( action != "do" ) return "";
  return cmd(i, arguments);
}

string CommandBase::type() const {
  return "Cm";
}

string CommandBase::doxygenType() const {
  return "Command";
}

CmdExUnknown::CmdExUnknown(const InterfaceBase & i,
			   const InterfacedBase & o, string c) {
  theMessage << "Could execute the command \"" << i.name()
	     << "\" for the object \"" << o.name() << "\" with argument " << c
	     << "because the command function threw an unknown exception.";
  severity(warning);
}

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "Command.tcc"
#endif

