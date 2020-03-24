// -*- C++ -*-
//
// Switch.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SwitchBase class.
//

#include "InterfacedBase.h"
#include "Switch.h"
#include "ThePEG/Utilities/StringUtils.h"

namespace ThePEG {

string SwitchBase::opttag(long opt) const {
  ostringstream ret;
  ret << opt;
  OptionMap::const_iterator oit = theOptions.find(opt);
  if ( oit ==  theOptions.end() )
    ret << " [Not a registered option] ";
  else
    ret << " [" << (oit->second).name() << "]";
  return ret.str();
}

string SwitchBase::exec(InterfacedBase & i, string action,
		    string arguments) const {
  ostringstream ret;

  if ( action == "get" ) {
    return opttag(get(i));
  }
  else if ( action == "def" ) {
    return opttag(def(i));
  }
  else if ( action == "setdef" &&
	    objectDefaults(i).find(name()) == objectDefaults(i).end() ) {
    setDef(i);
  } else if ( action == "set" || action == "newdef"  || action == "setdef" ) {
    if ( action == "setdef" ) arguments = objectDefaults(i)[name()];
    istringstream arg(arguments);
    long val;
    arg >> val;
    if ( !arg || theOptions.find(val) == theOptions.end() ) {
      string sval = StringUtils::car(arguments);
      StringMap::const_iterator sit = theOptionNames.find(sval);
      if ( sit == theOptionNames.end() ) {
	if ( sval == "true" && theOptions.find(1) != theOptions.end() )
	  val = 1;
	else if ( sval == "false" && theOptions.find(0) != theOptions.end() )
	  val = 0;
	else {
	  string errortext = "Error: no option '" + StringUtils::car(arguments) +
	    "' found for switch\n'" + i.fullName() + ":" +name() + "'\nValid options are: ";
	  for (StringMap::const_iterator it = theOptionNames.begin();
	       it != theOptionNames.end(); ++it )
	    errortext += it->first + ' ';
	  return errortext + '\n';
	}
      } else
	val = sit->second.value();
    }
    try {
      set(i, val);
      if ( action == "newdef" ) objectDefaults(i)[name()] = opttag(val);
      return ret.str();
    }
    catch ( Exception & ex ) {
      ex.handle();
    }
    istringstream arg2(arguments.c_str());
    string sval;
    arg2 >> sval;
    StringMap::const_iterator optit = theOptionNames.find(sval);
    if ( optit != theOptionNames.end() ) val = optit->second.value();
    set(i, val);
    if ( action == "newdef" ) objectDefaults(i)[name()] = opttag(val);
  }
  else if ( action == "notdef" ) {
    string deflt = opttag(def(i));
    if ( objectDefaults(i).find(name()) != objectDefaults(i).end() )
      deflt = objectDefaults(i)[name()];
    else if ( !hasDefault ) return "";
    if ( deflt != opttag(get(i)) ) return opttag(get(i)) + " (" + deflt + ")";
  }
  else
    throw InterExUnknown(*this, i);
  return ret.str();
}

string SwitchBase::type() const {
  return "Sw";
}

string SwitchBase::doxygenType() const {
  return "Switch";
}

string SwitchBase::fullDescription(const InterfacedBase & ib) const {
  ostringstream os;
  os << InterfaceBase::fullDescription(ib) << get(ib) << '\n'
     << def(ib) << '\n' << options().size() << '\n';
  for ( OptionMap::const_iterator it = options().begin();
	it != options().end(); ++it )
    os << it->second.value() << '\n'
       << it->second.name() << '\n'
       << it->second.description() << endl;
  return os .str();
}

SwExSetOpt::SwExSetOpt(const InterfaceBase & i,
		       const InterfacedBase & o, long v) {
  theMessage << "Could not set the switch \"" << i.name()
	     << "\" for the object \"" << o.name() << "\" to " << v
	     << "because it is not a registered option.";
  severity(setuperror);
}

SwExSetUnknown::SwExSetUnknown(const InterfaceBase & i,
			       const InterfacedBase & o, long v) {
  theMessage << "Could not set the switch \"" << i.name()
	     << "\" for the object \"" << o.name() << "\" to " << v
	     << " because the set function threw an unknown exception.";
  severity(setuperror);
}

SwExGetUnknown::SwExGetUnknown(const InterfaceBase & i,
			       const InterfacedBase & o, const char * s) {
  theMessage << "Could not get the " << s << " value of switch \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the get function threw an unknown exception.";
  severity(setuperror);
}

}

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "Switch.tcc"
#endif
