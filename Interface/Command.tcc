// -*- C++ -*-
//
// Command.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the Command class.
//

namespace ThePEG {

template <class T>
string Command<T>::cmd(InterfacedBase & i, string arg) const
  {
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theExeFn ) {
    try {
      string r = (t->*theExeFn)(arg);
      if ( r != "" ) i.touch();
      return r;
    }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw CmdExUnknown(*this, i, arg); }
  } else throw InterExSetup(*this, i);
}

}
