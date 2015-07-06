#ifndef fpudebug_H
#define fpudebug_H

#include <fpu_control.h>

/** Helper function to enable trapping of floating point exceptions. */
void fpudebug() {
    fpu_control_t cw;
  _FPU_GETCW(cw);
  cw &= ~(_FPU_MASK_OM);
  cw &= ~(_FPU_MASK_ZM);
  cw &= ~(_FPU_MASK_DM);
  cw &= ~(_FPU_MASK_IM);
  _FPU_SETCW(cw);
}

#endif
