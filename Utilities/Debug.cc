// -*- C++ -*-
//
// Debug.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Debug class.
//

#include "Debug.h"
#include "config.h"
#ifdef ThePEG_HAS_FPU_CONTROL
#include <fpu_control.h>
#endif
#ifdef ThePEG_HAS_FENV
#include <fenv.h>
#endif

using namespace ThePEG;

int Debug::level = 0;

bool Debug::isset = false;

std::vector<bool> Debug::debugItems;

void Debug::debugItem(int item, bool on) {
  if ( item < 0 ) return;
  debugItems.resize(item + 1, false);
  debugItems[item] = on;
}

void Debug::setDebug(int ilev) {
  if ( ilev < 0 ) debugItem(-ilev, true);
  else {
    level = ilev;
    isset = true;
  }
}

void Debug::unmaskFpuErrors() {
  unmaskFpuOverflow();
  unmaskFpuDivZero();
  //  unmaskFpuDenorm();
  unmaskFpuInvalid();
}

void Debug::unmaskFpuOverflow() {
#ifdef ThePEG_HAS_FENV
  feenableexcept(FE_OVERFLOW);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw &= ~(_FPU_MASK_OM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::unmaskFpuUnderflow() {
#ifdef ThePEG_HAS_FENV
  feenableexcept(FE_UNDERFLOW);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw &= ~(_FPU_MASK_UM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::unmaskFpuDivZero() {
#ifdef ThePEG_HAS_FENV
  feenableexcept(FE_DIVBYZERO);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw &= ~(_FPU_MASK_ZM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::unmaskFpuDenorm() {
#ifdef ThePEG_HAS_FENV
  feenableexcept(FE_INEXACT);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw &= ~(_FPU_MASK_DM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::unmaskFpuInvalid()  {
#ifdef ThePEG_HAS_FENV
  feenableexcept(FE_INVALID);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw &= ~(_FPU_MASK_IM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::maskFpuErrors() {
  maskFpuOverflow();
  maskFpuDivZero();
  maskFpuDenorm();
  maskFpuInvalid();
}

void Debug::maskFpuOverflow() {
#ifdef ThePEG_HAS_FENV
  fedisableexcept(FE_OVERFLOW);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw |= (_FPU_MASK_OM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::maskFpuUnderflow() {
#ifdef ThePEG_HAS_FENV
  fedisableexcept(FE_UNDERFLOW);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw |= (_FPU_MASK_UM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::maskFpuDivZero() {
#ifdef ThePEG_HAS_FENV
  fedisableexcept(FE_DIVBYZERO);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw |= (_FPU_MASK_ZM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::maskFpuDenorm() {
#ifdef ThePEG_HAS_FENV
  fedisableexcept(FE_INEXACT);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw |= (_FPU_MASK_DM);
  _FPU_SETCW(cw);
#endif
#endif
}

void Debug::maskFpuInvalid() {
#ifdef ThePEG_HAS_FENV
  fedisableexcept(FE_INVALID);
#else
#ifdef ThePEG_HAS_FPU_CONTROL
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw |= (_FPU_MASK_IM);
  _FPU_SETCW(cw);
#endif
#endif
}

