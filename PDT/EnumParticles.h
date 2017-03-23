// -*- C++ -*-
//
// EnumParticles.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_EnumParticles_H
#define ThePEG_EnumParticles_H

namespace ThePEG {


/**
 * The ParticleID namespace defines the ParticleCodes enumeration.
 */
namespace ParticleID {

/**
 * Enumeration to give identifiers to PDG id numbers.
 */
enum ParticleCodes {
  d = 1,
    /**< \f$\mbox{d}\f$
         (<code>d</code>) */
  dbar = -1,
    /**< \f$\overline{\mbox{d}}\f$
         (<code>dbar</code>) */
  u = 2,
    /**< \f$\mbox{u}\f$
         (<code>u</code>) */
  ubar = -2,
    /**< \f$\overline{\mbox{u}}\f$
         (<code>ubar</code>) */
  s = 3,
    /**< \f$\mbox{s}\f$
         (<code>s</code>) */
  sbar = -3,
    /**< \f$\overline{\mbox{s}}\f$
         (<code>sbar</code>) */
  c = 4,
    /**< \f$\mbox{c}\f$
         (<code>c</code>) */
  cbar = -4,
    /**< \f$\overline{\mbox{c}}\f$
         (<code>cbar</code>) */
  b = 5,
    /**< \f$\mbox{b}\f$
         (<code>b</code>) */
  bbar = -5,
    /**< \f$\overline{\mbox{b}}\f$
         (<code>bbar</code>) */
  t = 6,
    /**< \f$\mbox{t}\f$
         (<code>t</code>) */
  tbar = -6,
    /**< \f$\overline{\mbox{t}}\f$
         (<code>tbar</code>) */
  bprime = 7,
    /**< \f$\mbox{b}^{\prime }\f$
         (<code>b'</code>) */
  bprimebar = -7,
    /**< \f$\overline{\mbox{b}}^{\prime }\f$
         (<code>b'bar</code>) */
  tprime = 8,
    /**< \f$\mbox{t}^{\prime }\f$
         (<code>t'</code>) */
  tprimebar = -8,
    /**< \f$\overline{\mbox{t}}^{\prime }\f$
         (<code>t'bar</code>) */
  eminus = 11,
    /**< \f$\mbox{e}^{-}\f$
         (<code>e-</code>) */
  eplus = -11,
    /**< \f$\mbox{e}^{+}\f$
         (<code>e+</code>) */
  nu_e = 12,
    /**< \f$\nu _{e}\f$
         (<code>nu_e</code>) */
  nu_ebar = -12,
    /**< \f$\overline{\nu }_{e}\f$
         (<code>nu_ebar</code>) */
  muminus = 13,
    /**< \f$\mu ^{-}\f$
         (<code>mu-</code>) */
  muplus = -13,
    /**< \f$\mu ^{+}\f$
         (<code>mu+</code>) */
  nu_mu = 14,
    /**< \f$\nu _{\mu }\f$
         (<code>nu_mu</code>) */
  nu_mubar = -14,
    /**< \f$\overline{\nu }_{\mu }\f$
         (<code>nu_mubar</code>) */
  tauminus = 15,
    /**< \f$\tau ^{-}\f$
         (<code>tau-</code>) */
  tauplus = -15,
    /**< \f$\tau ^{+}\f$
         (<code>tau+</code>) */
  nu_tau = 16,
    /**< \f$\nu _{\tau }\f$
         (<code>nu_tau</code>) */
  nu_taubar = -16,
    /**< \f$\overline{\nu }_{\tau }\f$
         (<code>nu_taubar</code>) */
  tauprimeminus = 17,
    /**< \f$\tau ^{\prime -}\f$
         (<code>tau'-</code>) */
  tauprimeplus = -17,
    /**< \f$\tau ^{\prime +}\f$
         (<code>tau'+</code>) */
  nuprime_tau = 18,
    /**< \f$\nu ^{\prime }_{\tau }\f$
         (<code>nu'_tau</code>) */
  nuprime_taubar = -18,
    /**< \f$\overline{\nu }^{\prime }_{\tau }\f$
         (<code>nu'_taubar</code>) */
  g = 21,
    /**< \f$\mbox{g}\f$
         (<code>g</code>) */
  gamma = 22,
    /**< \f$\gamma \f$
         (<code>gamma</code>) */
  Z0 = 23,
    /**< \f$\mbox{Z}^{0 }\f$
         (<code>Z0</code>) */
  Wplus = 24,
    /**< \f$\mbox{W}^{+}\f$
         (<code>W+</code>) */
  Wminus = -24,
    /**< \f$\mbox{W}^{-}\f$
         (<code>W-</code>) */
  h0 = 25,
    /**< \f$\mbox{h}^{0 }\f$
         (<code>h0</code>) */
  Zprime0 = 32,
    /**< \f$\mbox{Z}^{\prime 0 }\f$
         (<code>Z'0</code>) */
  Zbis0 = 33,
    /**< \f$\mbox{Z}^{\prime\prime 0 }\f$
         (<code>Z"0</code>) */
  Wprimeplus = 34,
    /**< \f$\mbox{W}^{\prime +}\f$
         (<code>W'+</code>) */
  Wprimeminus = -34,
    /**< \f$\mbox{W}^{\prime -}\f$
         (<code>W'-</code>) */
  H0 = 35,
    /**< \f$\mbox{H}^{0 }\f$
         (<code>H0</code>) */
  A0 = 36,
    /**< \f$\mbox{A}^{0 }\f$
         (<code>A0</code>) */
  Hplus = 37,
    /**< \f$\mbox{H}^{+}\f$
         (<code>H+</code>) */
  Hminus = -37,
    /**< \f$\mbox{H}^{-}\f$
         (<code>H-</code>) */
  Graviton = 39,
    /**< \f${\cal G}\f$
         (<code>Graviton</code>) */
  R0 = 41,
    /**< \f$\mbox{R}^{0 }\f$
         (<code>R0</code>) */
  Rbar0 = -41,
    /**< \f$\overline{\mbox{R}}^{0 }\f$
         (<code>Rbar0</code>) */
  LQ_ue = 42,
    /**< \f$\mbox{L}_{Que}\f$
         (<code>LQ_ue</code>) */
  LQ_uebar = -42,
    /**< \f$\overline{\mbox{L}}_{Que}\f$
         (<code>LQ_uebar</code>) */
  reggeon = 110,
    /**< \f$I\!\!R\f$
         (<code>reggeon</code>) */
  pi0 = 111,
    /**< \f$\pi ^{0 }\f$
         (<code>pi0</code>) */
  rho0 = 113,
    /**< \f$\rho ^{0 }\f$
         (<code>rho0</code>) */
  a_20 = 115,
    /**< \f$\mbox{a}^{0 }_{2}\f$
         (<code>a_20</code>) */
  K_L0 = 130,
    /**< \f$\mbox{K}^{0 }_{L}\f$
         (<code>K_L0</code>) */
  piplus = 211,
    /**< \f$\pi ^{+}\f$
         (<code>pi+</code>) */
  piminus = -211,
    /**< \f$\pi ^{-}\f$
         (<code>pi-</code>) */
  rhoplus = 213,
    /**< \f$\rho ^{+}\f$
         (<code>rho+</code>) */
  rhominus = -213,
    /**< \f$\rho ^{-}\f$
         (<code>rho-</code>) */
  a_2plus = 215,
    /**< \f$\mbox{a}^{+}_{2}\f$
         (<code>a_2+</code>) */
  a_2minus = -215,
    /**< \f$\mbox{a}^{-}_{2}\f$
         (<code>a_2-</code>) */
  eta = 221,
    /**< \f$\eta \f$
         (<code>eta</code>) */
  omega = 223,
    /**< \f$\omega \f$
         (<code>omega</code>) */
  f_2 = 225,
    /**< \f$\mbox{f}_{2}\f$
         (<code>f_2</code>) */
  K_S0 = 310,
    /**< \f$\mbox{K}^{0 }_{S}\f$
         (<code>K_S0</code>) */
  K0 = 311,
    /**< \f$\mbox{K}^{0 }\f$
         (<code>K0</code>) */
  Kbar0 = -311,
    /**< \f$\overline{\mbox{K}}^{0 }\f$
         (<code>Kbar0</code>) */
  Kstar0 = 313,
    /**< \f$\mbox{K}^{*0 }\f$
         (<code>K*0</code>) */
  Kstarbar0 = -313,
    /**< \f$\overline{\mbox{K}}^{*0 }\f$
         (<code>K*bar0</code>) */
  Kstar_20 = 315,
    /**< \f$\mbox{K}^{*0 }_{2}\f$
         (<code>K*_20</code>) */
  Kstar_2bar0 = -315,
    /**< \f$\overline{\mbox{K}}^{*0 }_{2}\f$
         (<code>K*_2bar0</code>) */
  Kplus = 321,
    /**< \f$\mbox{K}^{+}\f$
         (<code>K+</code>) */
  Kminus = -321,
    /**< \f$\mbox{K}^{-}\f$
         (<code>K-</code>) */
  Kstarplus = 323,
    /**< \f$\mbox{K}^{*+}\f$
         (<code>K*+</code>) */
  Kstarminus = -323,
    /**< \f$\mbox{K}^{*-}\f$
         (<code>K*-</code>) */
  Kstar_2plus = 325,
    /**< \f$\mbox{K}^{*+}_{2}\f$
         (<code>K*_2+</code>) */
  Kstar_2minus = -325,
    /**< \f$\mbox{K}^{*-}_{2}\f$
         (<code>K*_2-</code>) */
  etaprime = 331,
    /**< \f$\eta ^{\prime }\f$
         (<code>eta'</code>) */
  phi = 333,
    /**< \f$\phi \f$
         (<code>phi</code>) */
  fprime_2 = 335,
    /**< \f$\mbox{f}^{\prime }_{2}\f$
         (<code>f'_2</code>) */
  Dplus = 411,
    /**< \f$\mbox{D}^{+}\f$
         (<code>D+</code>) */
  Dminus = -411,
    /**< \f$\mbox{D}^{-}\f$
         (<code>D-</code>) */
  Dstarplus = 413,
    /**< \f$\mbox{D}^{*+}\f$
         (<code>D*+</code>) */
  Dstarminus = -413,
    /**< \f$\mbox{D}^{*-}\f$
         (<code>D*-</code>) */
  Dstar_2plus = 415,
    /**< \f$\mbox{D}^{*+}_{2}\f$
         (<code>D*_2+</code>) */
  Dstar_2minus = -415,
    /**< \f$\mbox{D}^{*-}_{2}\f$
         (<code>D*_2-</code>) */
  D0 = 421,
    /**< \f$\mbox{D}^{0 }\f$
         (<code>D0</code>) */
  Dbar0 = -421,
    /**< \f$\overline{\mbox{D}}^{0 }\f$
         (<code>Dbar0</code>) */
  Dstar0 = 423,
    /**< \f$\mbox{D}^{*0 }\f$
         (<code>D*0</code>) */
  Dstarbar0 = -423,
    /**< \f$\overline{\mbox{D}}^{*0 }\f$
         (<code>D*bar0</code>) */
  Dstar_20 = 425,
    /**< \f$\mbox{D}^{*0 }_{2}\f$
         (<code>D*_20</code>) */
  Dstar_2bar0 = -425,
    /**< \f$\overline{\mbox{D}}^{*0 }_{2}\f$
         (<code>D*_2bar0</code>) */
  D_splus = 431,
    /**< \f$\mbox{D}^{+}_{s}\f$
         (<code>D_s+</code>) */
  D_sminus = -431,
    /**< \f$\mbox{D}^{-}_{s}\f$
         (<code>D_s-</code>) */
  Dstar_splus = 433,
    /**< \f$\mbox{D}^{*+}_{s}\f$
         (<code>D*_s+</code>) */
  Dstar_sminus = -433,
    /**< \f$\mbox{D}^{*-}_{s}\f$
         (<code>D*_s-</code>) */
  Dstar_2splus = 435,
    /**< \f$\mbox{D}^{*+}_{2s}\f$
         (<code>D*_2s+</code>) */
  Dstar_2sminus = -435,
    /**< \f$\mbox{D}^{*-}_{2s}\f$
         (<code>D*_2s-</code>) */
  eta_c = 441,
    /**< \f$\eta _{c}\f$
         (<code>eta_c</code>) */
  Jpsi = 443,
    /**< \f$J/\psi \f$
         (<code>J/psi</code>) */
  chi_2c = 445,
    /**< \f$\chi _{2c}\f$
         (<code>chi_2c</code>) */
  B0 = 511,
    /**< \f$\mbox{B}^{0 }\f$
         (<code>B0</code>) */
  Bbar0 = -511,
    /**< \f$\overline{\mbox{B}}^{0 }\f$
         (<code>Bbar0</code>) */
  Bstar0 = 513,
    /**< \f$\mbox{B}^{*0 }\f$
         (<code>B*0</code>) */
  Bstarbar0 = -513,
    /**< \f$\overline{\mbox{B}}^{*0 }\f$
         (<code>B*bar0</code>) */
  Bstar_20 = 515,
    /**< \f$\mbox{B}^{*0 }_{2}\f$
         (<code>B*_20</code>) */
  Bstar_2bar0 = -515,
    /**< \f$\overline{\mbox{B}}^{*0 }_{2}\f$
         (<code>B*_2bar0</code>) */
  Bplus = 521,
    /**< \f$\mbox{B}^{+}\f$
         (<code>B+</code>) */
  Bminus = -521,
    /**< \f$\mbox{B}^{-}\f$
         (<code>B-</code>) */
  Bstarplus = 523,
    /**< \f$\mbox{B}^{*+}\f$
         (<code>B*+</code>) */
  Bstarminus = -523,
    /**< \f$\mbox{B}^{*-}\f$
         (<code>B*-</code>) */
  Bstar_2plus = 525,
    /**< \f$\mbox{B}^{*+}_{2}\f$
         (<code>B*_2+</code>) */
  Bstar_2minus = -525,
    /**< \f$\mbox{B}^{*-}_{2}\f$
         (<code>B*_2-</code>) */
  B_s0 = 531,
    /**< \f$\mbox{B}^{0 }_{s}\f$
         (<code>B_s0</code>) */
  B_sbar0 = -531,
    /**< \f$\overline{\mbox{B}}^{0 }_{s}\f$
         (<code>B_sbar0</code>) */
  Bstar_s0 = 533,
    /**< \f$\mbox{B}^{*0 }_{s}\f$
         (<code>B*_s0</code>) */
  Bstar_sbar0 = -533,
    /**< \f$\overline{\mbox{B}}^{*0 }_{s}\f$
         (<code>B*_sbar0</code>) */
  Bstar_2s0 = 535,
    /**< \f$\mbox{B}^{*0 }_{2s}\f$
         (<code>B*_2s0</code>) */
  Bstar_2sbar0 = -535,
    /**< \f$\overline{\mbox{B}}^{*0 }_{2s}\f$
         (<code>B*_2sbar0</code>) */
  B_cplus = 541,
    /**< \f$\mbox{B}^{+}_{c}\f$
         (<code>B_c+</code>) */
  B_cminus = -541,
    /**< \f$\mbox{B}^{-}_{c}\f$
         (<code>B_c-</code>) */
  Bstar_cplus = 543,
    /**< \f$\mbox{B}^{*+}_{c}\f$
         (<code>B*_c+</code>) */
  Bstar_cminus = -543,
    /**< \f$\mbox{B}^{*-}_{c}\f$
         (<code>B*_c-</code>) */
  Bstar_2cplus = 545,
    /**< \f$\mbox{B}^{*+}_{2c}\f$
         (<code>B*_2c+</code>) */
  Bstar_2cminus = -545,
    /**< \f$\mbox{B}^{*-}_{2c}\f$
         (<code>B*_2c-</code>) */
  eta_b = 551,
    /**< \f$\eta _{b}\f$
         (<code>eta_b</code>) */
  Upsilon = 553,
    /**< \f$\Upsilon  \f$
         (<code>Upsilon</code>) */
  chi_2b = 555,
    /**< \f$\chi _{2b}\f$
         (<code>chi_2b</code>) */
  pomeron = 990,
    /**< \f$I\!\!P\f$
         (<code>pomeron</code>) */
  dd_1 = 1103,
    /**< \f$\mbox{dd}_{1}\f$
         (<code>dd_1</code>) */
  dd_1bar = -1103,
    /**< \f$\overline{\mbox{dd}}_{1}\f$
         (<code>dd_1bar</code>) */
  Deltaminus = 1114,
    /**< \f$\Delta ^{-}\f$
         (<code>Delta-</code>) */
  Deltabarplus = -1114,
    /**< \f$\overline{\Delta }^{+}\f$
         (<code>Deltabar+</code>) */
  ud_0 = 2101,
    /**< \f$\mbox{ud}^{0 }\f$
         (<code>ud_0</code>) */
  ud_0bar = -2101,
    /**< \f$\overline{\mbox{ud}}^{0 }\f$
         (<code>ud_0bar</code>) */
  ud_1 = 2103,
    /**< \f$\mbox{ud}_{1}\f$
         (<code>ud_1</code>) */
  ud_1bar = -2103,
    /**< \f$\overline{\mbox{ud}}_{1}\f$
         (<code>ud_1bar</code>) */
  n0 = 2112,
    /**< \f$\mbox{n}^{0 }\f$
         (<code>n0</code>) */
  nbar0 = -2112,
    /**< \f$\overline{\mbox{n}}^{0 }\f$
         (<code>nbar0</code>) */
  Delta0 = 2114,
    /**< \f$\Delta ^{0 }\f$
         (<code>Delta0</code>) */
  Deltabar0 = -2114,
    /**< \f$\overline{\Delta }^{0 }\f$
         (<code>Deltabar0</code>) */
  uu_1 = 2203,
    /**< \f$\mbox{uu}_{1}\f$
         (<code>uu_1</code>) */
  uu_1bar = -2203,
    /**< \f$\overline{\mbox{uu}}_{1}\f$
         (<code>uu_1bar</code>) */
  pplus = 2212,
    /**< \f$\mbox{p}^{+}\f$
         (<code>p+</code>) */
  pbarminus = -2212,
    /**< \f$\overline{\mbox{p}}^{-}\f$
         (<code>pbar-</code>) */
  Deltaplus = 2214,
    /**< \f$\Delta ^{+}\f$
         (<code>Delta+</code>) */
  Deltabarminus = -2214,
    /**< \f$\overline{\Delta }^{-}\f$
         (<code>Deltabar-</code>) */
  Deltaplus2 = 2224,
    /**< \f$\Delta ^{++}\f$
         (<code>Delta++</code>) */
  Deltabarminus2 = -2224,
    /**< \f$\overline{\Delta }^{--}\f$
         (<code>Deltabar--</code>) */
  sd_0 = 3101,
    /**< \f$\mbox{sd}^{0 }\f$
         (<code>sd_0</code>) */
  sd_0bar = -3101,
    /**< \f$\overline{\mbox{sd}}^{0 }\f$
         (<code>sd_0bar</code>) */
  sd_1 = 3103,
    /**< \f$\mbox{sd}_{1}\f$
         (<code>sd_1</code>) */
  sd_1bar = -3103,
    /**< \f$\overline{\mbox{sd}}_{1}\f$
         (<code>sd_1bar</code>) */
  Sigmaminus = 3112,
    /**< \f$\Sigma ^{-}\f$
         (<code>Sigma-</code>) */
  Sigmabarplus = -3112,
    /**< \f$\overline{\Sigma }^{+}\f$
         (<code>Sigmabar+</code>) */
  Sigmastarminus = 3114,
    /**< \f$\Sigma ^{*-}\f$
         (<code>Sigma*-</code>) */
  Sigmastarbarplus = -3114,
    /**< \f$\overline{\Sigma }^{*+}\f$
         (<code>Sigma*bar+</code>) */
  Lambda0 = 3122,
    /**< \f$\Lambda ^{0 }\f$
         (<code>Lambda0</code>) */
  Lambdabar0 = -3122,
    /**< \f$\overline{\Lambda }^{0 }\f$
         (<code>Lambdabar0</code>) */
  su_0 = 3201,
    /**< \f$\mbox{su}^{0 }\f$
         (<code>su_0</code>) */
  su_0bar = -3201,
    /**< \f$\overline{\mbox{su}}^{0 }\f$
         (<code>su_0bar</code>) */
  su_1 = 3203,
    /**< \f$\mbox{su}_{1}\f$
         (<code>su_1</code>) */
  su_1bar = -3203,
    /**< \f$\overline{\mbox{su}}_{1}\f$
         (<code>su_1bar</code>) */
  Sigma0 = 3212,
    /**< \f$\Sigma ^{0 }\f$
         (<code>Sigma0</code>) */
  Sigmabar0 = -3212,
    /**< \f$\overline{\Sigma }^{0 }\f$
         (<code>Sigmabar0</code>) */
  Sigmastar0 = 3214,
    /**< \f$\Sigma ^{*0 }\f$
         (<code>Sigma*0</code>) */
  Sigmastarbar0 = -3214,
    /**< \f$\overline{\Sigma }^{*0 }\f$
         (<code>Sigma*bar0</code>) */
  Sigmaplus = 3222,
    /**< \f$\Sigma ^{+}\f$
         (<code>Sigma+</code>) */
  Sigmabarminus = -3222,
    /**< \f$\overline{\Sigma }^{-}\f$
         (<code>Sigmabar-</code>) */
  Sigmastarplus = 3224,
    /**< \f$\Sigma ^{*+}\f$
         (<code>Sigma*+</code>) */
  Sigmastarbarminus = -3224,
    /**< \f$\overline{\Sigma }^{*-}\f$
         (<code>Sigma*bar-</code>) */
  ss_1 = 3303,
    /**< \f$\mbox{ss}_{1}\f$
         (<code>ss_1</code>) */
  ss_1bar = -3303,
    /**< \f$\overline{\mbox{ss}}_{1}\f$
         (<code>ss_1bar</code>) */
  Ximinus = 3312,
    /**< \f$\Xi ^{-}\f$
         (<code>Xi-</code>) */
  Xibarplus = -3312,
    /**< \f$\overline{\Xi }^{+}\f$
         (<code>Xibar+</code>) */
  Xistarminus = 3314,
    /**< \f$\Xi ^{*-}\f$
         (<code>Xi*-</code>) */
  Xistarbarplus = -3314,
    /**< \f$\overline{\Xi }^{*+}\f$
         (<code>Xi*bar+</code>) */
  Xi0 = 3322,
    /**< \f$\Xi ^{0 }\f$
         (<code>Xi0</code>) */
  Xibar0 = -3322,
    /**< \f$\overline{\Xi }^{0 }\f$
         (<code>Xibar0</code>) */
  Xistar0 = 3324,
    /**< \f$\Xi ^{*0 }\f$
         (<code>Xi*0</code>) */
  Xistarbar0 = -3324,
    /**< \f$\overline{\Xi }^{*0 }\f$
         (<code>Xi*bar0</code>) */
  Omegaminus = 3334,
    /**< \f$\Omega ^{-}\f$
         (<code>Omega-</code>) */
  Omegabarplus = -3334,
    /**< \f$\overline{\Omega }^{+}\f$
         (<code>Omegabar+</code>) */
  cd_0 = 4101,
    /**< \f$\mbox{cd}^{0 }\f$
         (<code>cd_0</code>) */
  cd_0bar = -4101,
    /**< \f$\overline{\mbox{cd}}^{0 }\f$
         (<code>cd_0bar</code>) */
  cd_1 = 4103,
    /**< \f$\mbox{cd}_{1}\f$
         (<code>cd_1</code>) */
  cd_1bar = -4103,
    /**< \f$\overline{\mbox{cd}}_{1}\f$
         (<code>cd_1bar</code>) */
  Sigma_c0 = 4112,
    /**< \f$\Sigma ^{0 }_{c}\f$
         (<code>Sigma_c0</code>) */
  Sigma_cbar0 = -4112,
    /**< \f$\overline{\Sigma }^{0 }_{c}\f$
         (<code>Sigma_cbar0</code>) */
  Sigmastar_c0 = 4114,
    /**< \f$\Sigma ^{*0 }_{c}\f$
         (<code>Sigma*_c0</code>) */
  Sigmastar_cbar0 = -4114,
    /**< \f$\overline{\Sigma }^{*0 }_{c}\f$
         (<code>Sigma*_cbar0</code>) */
  Lambda_cplus = 4122,
    /**< \f$\Lambda ^{+}_{c}\f$
         (<code>Lambda_c+</code>) */
  Lambda_cbarminus = -4122,
    /**< \f$\overline{\Lambda }^{-}_{c}\f$
         (<code>Lambda_cbar-</code>) */
  Xi_c0 = 4132,
    /**< \f$\Xi ^{0 }_{c}\f$
         (<code>Xi_c0</code>) */
  Xi_cbar0 = -4132,
    /**< \f$\overline{\Xi }^{0 }_{c}\f$
         (<code>Xi_cbar0</code>) */
  cu_0 = 4201,
    /**< \f$\mbox{cu}^{0 }\f$
         (<code>cu_0</code>) */
  cu_0bar = -4201,
    /**< \f$\overline{\mbox{cu}}^{0 }\f$
         (<code>cu_0bar</code>) */
  cu_1 = 4203,
    /**< \f$\mbox{cu}_{1}\f$
         (<code>cu_1</code>) */
  cu_1bar = -4203,
    /**< \f$\overline{\mbox{cu}}_{1}\f$
         (<code>cu_1bar</code>) */
  Sigma_cplus = 4212,
    /**< \f$\Sigma ^{+}_{c}\f$
         (<code>Sigma_c+</code>) */
  Sigma_cbarminus = -4212,
    /**< \f$\overline{\Sigma }^{-}_{c}\f$
         (<code>Sigma_cbar-</code>) */
  Sigmastar_cplus = 4214,
    /**< \f$\Sigma ^{*+}_{c}\f$
         (<code>Sigma*_c+</code>) */
  Sigmastar_cbarminus = -4214,
    /**< \f$\overline{\Sigma }^{*-}_{c}\f$
         (<code>Sigma*_cbar-</code>) */
  Sigma_cplus2 = 4222,
    /**< \f$\Sigma ^{++}_{c}\f$
         (<code>Sigma_c++</code>) */
  Sigma_cbarminus2 = -4222,
    /**< \f$\overline{\Sigma }^{--}_{c}\f$
         (<code>Sigma_cbar--</code>) */
  Sigmastar_cplus2 = 4224,
    /**< \f$\Sigma ^{*++}_{c}\f$
         (<code>Sigma*_c++</code>) */
  Sigmastar_cbarminus2 = -4224,
    /**< \f$\overline{\Sigma }^{*--}_{c}\f$
         (<code>Sigma*_cbar--</code>) */
  Xi_cplus = 4232,
    /**< \f$\Xi ^{+}_{c}\f$
         (<code>Xi_c+</code>) */
  Xi_cbarminus = -4232,
    /**< \f$\overline{\Xi }^{-}_{c}\f$
         (<code>Xi_cbar-</code>) */
  cs_0 = 4301,
    /**< \f$\mbox{cs}^{0 }\f$
         (<code>cs_0</code>) */
  cs_0bar = -4301,
    /**< \f$\overline{\mbox{cs}}^{0 }\f$
         (<code>cs_0bar</code>) */
  cs_1 = 4303,
    /**< \f$\mbox{cs}_{1}\f$
         (<code>cs_1</code>) */
  cs_1bar = -4303,
    /**< \f$\overline{\mbox{cs}}_{1}\f$
         (<code>cs_1bar</code>) */
  Xiprime_c0 = 4312,
    /**< \f$\Xi ^{\prime 0 }_{c}\f$
         (<code>Xi'_c0</code>) */
  Xiprime_cbar0 = -4312,
    /**< \f$\overline{\Xi }^{\prime 0 }_{c}\f$
         (<code>Xi'_cbar0</code>) */
  Xistar_c0 = 4314,
    /**< \f$\Xi ^{*0 }_{c}\f$
         (<code>Xi*_c0</code>) */
  Xistar_cbar0 = -4314,
    /**< \f$\overline{\Xi }^{*0 }_{c}\f$
         (<code>Xi*_cbar0</code>) */
  Xiprime_cplus = 4322,
    /**< \f$\Xi ^{\prime +}_{c}\f$
         (<code>Xi'_c+</code>) */
  Xiprime_cbarminus = -4322,
    /**< \f$\overline{\Xi }^{\prime -}_{c}\f$
         (<code>Xi'_cbar-</code>) */
  Xistar_cplus = 4324,
    /**< \f$\Xi ^{*+}_{c}\f$
         (<code>Xi*_c+</code>) */
  Xistar_cbarminus = -4324,
    /**< \f$\overline{\Xi }^{*-}_{c}\f$
         (<code>Xi*_cbar-</code>) */
  Omega_c0 = 4332,
    /**< \f$\Omega ^{0 }_{c}\f$
         (<code>Omega_c0</code>) */
  Omega_cbar0 = -4332,
    /**< \f$\overline{\Omega }^{0 }_{c}\f$
         (<code>Omega_cbar0</code>) */
  Omegastar_c0 = 4334,
    /**< \f$\Omega ^{*0 }_{c}\f$
         (<code>Omega*_c0</code>) */
  Omegastar_cbar0 = -4334,
    /**< \f$\overline{\Omega }^{*0 }_{c}\f$
         (<code>Omega*_cbar0</code>) */
  cc_1 = 4403,
    /**< \f$\mbox{cc}_{1}\f$
         (<code>cc_1</code>) */
  cc_1bar = -4403,
    /**< \f$\overline{\mbox{cc}}_{1}\f$
         (<code>cc_1bar</code>) */
  Xi_ccplus = 4412,
    /**< \f$\Xi ^{+}_{cc}\f$
         (<code>Xi_cc+</code>) */
  Xi_ccbarminus = -4412,
    /**< \f$\overline{\Xi }^{-}_{cc}\f$
         (<code>Xi_ccbar-</code>) */
  Xistar_ccplus = 4414,
    /**< \f$\Xi ^{*+}_{cc}\f$
         (<code>Xi*_cc+</code>) */
  Xistar_ccbarminus = -4414,
    /**< \f$\overline{\Xi }^{*-}_{cc}\f$
         (<code>Xi*_ccbar-</code>) */
  Xi_ccplus2 = 4422,
    /**< \f$\Xi ^{++}_{cc}\f$
         (<code>Xi_cc++</code>) */
  Xi_ccbarminus2 = -4422,
    /**< \f$\overline{\Xi }^{--}_{cc}\f$
         (<code>Xi_ccbar--</code>) */
  Xistar_ccplus2 = 4424,
    /**< \f$\Xi ^{*++}_{cc}\f$
         (<code>Xi*_cc++</code>) */
  Xistar_ccbarminus2 = -4424,
    /**< \f$\overline{\Xi }^{*--}_{cc}\f$
         (<code>Xi*_ccbar--</code>) */
  Omega_ccplus = 4432,
    /**< \f$\Omega ^{+}_{cc}\f$
         (<code>Omega_cc+</code>) */
  Omega_ccbarminus = -4432,
    /**< \f$\overline{\Omega }^{-}_{cc}\f$
         (<code>Omega_ccbar-</code>) */
  Omegastar_ccplus = 4434,
    /**< \f$\Omega ^{*+}_{cc}\f$
         (<code>Omega*_cc+</code>) */
  Omegastar_ccbarminus = -4434,
    /**< \f$\overline{\Omega }^{*-}_{cc}\f$
         (<code>Omega*_ccbar-</code>) */
  Omegastar_cccplus2 = 4444,
    /**< \f$\Omega ^{*++}_{ccc}\f$
         (<code>Omega*_ccc++</code>) */
  Omegastar_cccbarminus = -4444,
    /**< \f$\overline{\Omega }^{*-}_{ccc}\f$
         (<code>Omega*_cccbar-</code>) */
  bd_0 = 5101,
    /**< \f$\mbox{bd}^{0 }\f$
         (<code>bd_0</code>) */
  bd_0bar = -5101,
    /**< \f$\overline{\mbox{bd}}^{0 }\f$
         (<code>bd_0bar</code>) */
  bd_1 = 5103,
    /**< \f$\mbox{bd}_{1}\f$
         (<code>bd_1</code>) */
  bd_1bar = -5103,
    /**< \f$\overline{\mbox{bd}}_{1}\f$
         (<code>bd_1bar</code>) */
  Sigma_bminus = 5112,
    /**< \f$\Sigma ^{-}_{b}\f$
         (<code>Sigma_b-</code>) */
  Sigma_bbarplus = -5112,
    /**< \f$\overline{\Sigma }^{+}_{b}\f$
         (<code>Sigma_bbar+</code>) */
  Sigmastar_bminus = 5114,
    /**< \f$\Sigma ^{*-}_{b}\f$
         (<code>Sigma*_b-</code>) */
  Sigmastar_bbarplus = -5114,
    /**< \f$\overline{\Sigma }^{*+}_{b}\f$
         (<code>Sigma*_bbar+</code>) */
  Lambda_b0 = 5122,
    /**< \f$\Lambda ^{0 }_{b}\f$
         (<code>Lambda_b0</code>) */
  Lambda_bbar0 = -5122,
    /**< \f$\overline{\Lambda }^{0 }_{b}\f$
         (<code>Lambda_bbar0</code>) */
  Xi_bminus = 5132,
    /**< \f$\Xi ^{-}_{b}\f$
         (<code>Xi_b-</code>) */
  Xi_bbarplus = -5132,
    /**< \f$\overline{\Xi }^{+}_{b}\f$
         (<code>Xi_bbar+</code>) */
  Xi_bc0 = 5142,
    /**< \f$\Xi ^{0 }_{bc}\f$
         (<code>Xi_bc0</code>) */
  Xi_bcbar0 = -5142,
    /**< \f$\overline{\Xi }^{0 }_{bc}\f$
         (<code>Xi_bcbar0</code>) */
  bu_0 = 5201,
    /**< \f$\mbox{bu}^{0 }\f$
         (<code>bu_0</code>) */
  bu_0bar = -5201,
    /**< \f$\overline{\mbox{bu}}^{0 }\f$
         (<code>bu_0bar</code>) */
  bu_1 = 5203,
    /**< \f$\mbox{bu}_{1}\f$
         (<code>bu_1</code>) */
  bu_1bar = -5203,
    /**< \f$\overline{\mbox{bu}}_{1}\f$
         (<code>bu_1bar</code>) */
  Sigma_b0 = 5212,
    /**< \f$\Sigma ^{0 }_{b}\f$
         (<code>Sigma_b0</code>) */
  Sigma_bbar0 = -5212,
    /**< \f$\overline{\Sigma }^{0 }_{b}\f$
         (<code>Sigma_bbar0</code>) */
  Sigmastar_b0 = 5214,
    /**< \f$\Sigma ^{*0 }_{b}\f$
         (<code>Sigma*_b0</code>) */
  Sigmastar_bbar0 = -5214,
    /**< \f$\overline{\Sigma }^{*0 }_{b}\f$
         (<code>Sigma*_bbar0</code>) */
  Sigma_bplus = 5222,
    /**< \f$\Sigma ^{+}_{b}\f$
         (<code>Sigma_b+</code>) */
  Sigma_bbarminus = -5222,
    /**< \f$\overline{\Sigma }^{-}_{b}\f$
         (<code>Sigma_bbar-</code>) */
  Sigmastar_bplus = 5224,
    /**< \f$\Sigma ^{*+}_{b}\f$
         (<code>Sigma*_b+</code>) */
  Sigmastar_bbarminus = -5224,
    /**< \f$\overline{\Sigma }^{*-}_{b}\f$
         (<code>Sigma*_bbar-</code>) */
  Xi_b0 = 5232,
    /**< \f$\Xi ^{0 }_{b}\f$
         (<code>Xi_b0</code>) */
  Xi_bbar0 = -5232,
    /**< \f$\overline{\Xi }^{0 }_{b}\f$
         (<code>Xi_bbar0</code>) */
  Xi_bcplus = 5242,
    /**< \f$\Xi ^{+}_{bc}\f$
         (<code>Xi_bc+</code>) */
  Xi_bcbarminus = -5242,
    /**< \f$\overline{\Xi }^{-}_{bc}\f$
         (<code>Xi_bcbar-</code>) */
  bs_0 = 5301,
    /**< \f$\mbox{bs}^{0 }\f$
         (<code>bs_0</code>) */
  bs_0bar = -5301,
    /**< \f$\overline{\mbox{bs}}^{0 }\f$
         (<code>bs_0bar</code>) */
  bs_1 = 5303,
    /**< \f$\mbox{bs}_{1}\f$
         (<code>bs_1</code>) */
  bs_1bar = -5303,
    /**< \f$\overline{\mbox{bs}}_{1}\f$
         (<code>bs_1bar</code>) */
  Xiprime_bminus = 5312,
    /**< \f$\Xi ^{\prime -}_{b}\f$
         (<code>Xi'_b-</code>) */
  Xiprime_bbarplus = -5312,
    /**< \f$\overline{\Xi }^{\prime +}_{b}\f$
         (<code>Xi'_bbar+</code>) */
  Xistar_bminus = 5314,
    /**< \f$\Xi ^{*-}_{b}\f$
         (<code>Xi*_b-</code>) */
  Xistar_bbarplus = -5314,
    /**< \f$\overline{\Xi }^{*+}_{b}\f$
         (<code>Xi*_bbar+</code>) */
  Xiprime_b0 = 5322,
    /**< \f$\Xi ^{\prime 0 }_{b}\f$
         (<code>Xi'_b0</code>) */
  Xiprime_bbar0 = -5322,
    /**< \f$\overline{\Xi }^{\prime 0 }_{b}\f$
         (<code>Xi'_bbar0</code>) */
  Xistar_b0 = 5324,
    /**< \f$\Xi ^{*0 }_{b}\f$
         (<code>Xi*_b0</code>) */
  Xistar_bbar0 = -5324,
    /**< \f$\overline{\Xi }^{*0 }_{b}\f$
         (<code>Xi*_bbar0</code>) */
  Omega_bminus = 5332,
    /**< \f$\Omega ^{-}_{b}\f$
         (<code>Omega_b-</code>) */
  Omega_bbarplus = -5332,
    /**< \f$\overline{\Omega }^{+}_{b}\f$
         (<code>Omega_bbar+</code>) */
  Omegastar_bminus = 5334,
    /**< \f$\Omega ^{*-}_{b}\f$
         (<code>Omega*_b-</code>) */
  Omegastar_bbarplus = -5334,
    /**< \f$\overline{\Omega }^{*+}_{b}\f$
         (<code>Omega*_bbar+</code>) */
  Omega_bc0 = 5342,
    /**< \f$\Omega ^{0 }_{bc}\f$
         (<code>Omega_bc0</code>) */
  Omega_bcbar0 = -5342,
    /**< \f$\overline{\Omega }^{0 }_{bc}\f$
         (<code>Omega_bcbar0</code>) */
  bc_0 = 5401,
    /**< \f$\mbox{bc}^{0 }\f$
         (<code>bc_0</code>) */
  bc_0bar = -5401,
    /**< \f$\overline{\mbox{bc}}^{0 }\f$
         (<code>bc_0bar</code>) */
  bc_1 = 5403,
    /**< \f$\mbox{bc}_{1}\f$
         (<code>bc_1</code>) */
  bc_1bar = -5403,
    /**< \f$\overline{\mbox{bc}}_{1}\f$
         (<code>bc_1bar</code>) */
  Xiprime_bc0 = 5412,
    /**< \f$\Xi ^{\prime 0 }_{bc}\f$
         (<code>Xi'_bc0</code>) */
  Xiprime_bcbar0 = -5412,
    /**< \f$\overline{\Xi }^{\prime 0 }_{bc}\f$
         (<code>Xi'_bcbar0</code>) */
  Xistar_bc0 = 5414,
    /**< \f$\Xi ^{*0 }_{bc}\f$
         (<code>Xi*_bc0</code>) */
  Xistar_bcbar0 = -5414,
    /**< \f$\overline{\Xi }^{*0 }_{bc}\f$
         (<code>Xi*_bcbar0</code>) */
  Xiprime_bcplus = 5422,
    /**< \f$\Xi ^{\prime +}_{bc}\f$
         (<code>Xi'_bc+</code>) */
  Xiprime_bcbarminus = -5422,
    /**< \f$\overline{\Xi }^{\prime -}_{bc}\f$
         (<code>Xi'_bcbar-</code>) */
  Xistar_bcplus = 5424,
    /**< \f$\Xi ^{*+}_{bc}\f$
         (<code>Xi*_bc+</code>) */
  Xistar_bcbarminus = -5424,
    /**< \f$\overline{\Xi }^{*-}_{bc}\f$
         (<code>Xi*_bcbar-</code>) */
  Omegaprime_bc0 = 5432,
    /**< \f$\Omega ^{\prime 0 }_{bc}\f$
         (<code>Omega'_bc0</code>) */
  Omegaprime_bcba = -5432,
    /**< \f$\overline{\Omega }^{\prime }_{bc}\f$
         (<code>Omega'_bcba</code>) */
  Omegastar_bc0 = 5434,
    /**< \f$\Omega ^{*0 }_{bc}\f$
         (<code>Omega*_bc0</code>) */
  Omegastar_bcbar0 = -5434,
    /**< \f$\overline{\Omega }^{*0 }_{bc}\f$
         (<code>Omega*_bcbar0</code>) */
  Omega_bccplus = 5442,
    /**< \f$\Omega ^{+}_{bcc}\f$
         (<code>Omega_bcc+</code>) */
  Omega_bccbarminus = -5442,
    /**< \f$\overline{\Omega }^{-}_{bcc}\f$
         (<code>Omega_bccbar-</code>) */
  Omegastar_bccplus = 5444,
    /**< \f$\Omega ^{*+}_{bcc}\f$
         (<code>Omega*_bcc+</code>) */
  Omegastar_bccbarminus = -5444,
    /**< \f$\overline{\Omega }^{*-}_{bcc}\f$
         (<code>Omega*_bccbar-</code>) */
  bb_1 = 5503,
    /**< \f$\mbox{bb}_{1}\f$
         (<code>bb_1</code>) */
  bb_1bar = -5503,
    /**< \f$\overline{\mbox{bb}}_{1}\f$
         (<code>bb_1bar</code>) */
  Xi_bbminus = 5512,
    /**< \f$\Xi ^{-}_{bb}\f$
         (<code>Xi_bb-</code>) */
  Xi_bbbarplus = -5512,
    /**< \f$\overline{\Xi }^{+}_{bb}\f$
         (<code>Xi_bbbar+</code>) */
  Xistar_bbminus = 5514,
    /**< \f$\Xi ^{*-}_{bb}\f$
         (<code>Xi*_bb-</code>) */
  Xistar_bbbarplus = -5514,
    /**< \f$\overline{\Xi }^{*+}_{bb}\f$
         (<code>Xi*_bbbar+</code>) */
  Xi_bb0 = 5522,
    /**< \f$\Xi ^{0 }_{bb}\f$
         (<code>Xi_bb0</code>) */
  Xi_bbbar0 = -5522,
    /**< \f$\overline{\Xi }^{0 }_{bb}\f$
         (<code>Xi_bbbar0</code>) */
  Xistar_bb0 = 5524,
    /**< \f$\Xi ^{*0 }_{bb}\f$
         (<code>Xi*_bb0</code>) */
  Xistar_bbbar0 = -5524,
    /**< \f$\overline{\Xi }^{*0 }_{bb}\f$
         (<code>Xi*_bbbar0</code>) */
  Omega_bbminus = 5532,
    /**< \f$\Omega ^{-}_{bb}\f$
         (<code>Omega_bb-</code>) */
  Omega_bbbarplus = -5532,
    /**< \f$\overline{\Omega }^{+}_{bb}\f$
         (<code>Omega_bbbar+</code>) */
  Omegastar_bbminus = 5534,
    /**< \f$\Omega ^{*-}_{bb}\f$
         (<code>Omega*_bb-</code>) */
  Omegastar_bbbarplus = -5534,
    /**< \f$\overline{\Omega }^{*+}_{bb}\f$
         (<code>Omega*_bbbar+</code>) */
  Omega_bbc0 = 5542,
    /**< \f$\Omega ^{0 }_{bbc}\f$
         (<code>Omega_bbc0</code>) */
  Omega_bbcbar0 = -5542,
    /**< \f$\overline{\Omega }^{0 }_{bbc}\f$
         (<code>Omega_bbcbar0</code>) */
  Omegastar_bbc0 = 5544,
    /**< \f$\Omega ^{*0 }_{bbc}\f$
         (<code>Omega*_bbc0</code>) */
  Omegastar_bbcbar0 = -5544,
    /**< \f$\overline{\Omega }^{*0 }_{bbc}\f$
         (<code>Omega*_bbcbar0</code>) */
  Omegastar_bbbminus = 5554,
    /**< \f$\Omega ^{*-}_{bbb}\f$
         (<code>Omega*_bbb-</code>) */
  Omegastar_bbbbarplus = -5554,
    /**< \f$\overline{\Omega }^{*+}_{bbb}\f$
         (<code>Omega*_bbbbar+</code>) */
  a_00 = 9000111,
    /**< \f$\mbox{a}^{0 }\f$
         (<code>a_00</code>) */
  b_10 = 10113,
    /**< \f$\mbox{b}^{0 }_{1}\f$
         (<code>b_10</code>) */
  a_0plus = 9000211,
    /**< \f$\mbox{a}^{0 +}\f$
         (<code>a_0+</code>) */
  a_0minus = -9000211,
    /**< \f$\mbox{a}^{0 -}\f$
         (<code>a_0-</code>) */
  b_1plus = 10213,
    /**< \f$\mbox{b}^{+}_{1}\f$
         (<code>b_1+</code>) */
  b_1minus = -10213,
    /**< \f$\mbox{b}^{-}_{1}\f$
         (<code>b_1-</code>) */
  f_0 = 9010221,
    /**< \f$\mbox{f}^{0 }\f$
         (<code>f_0</code>) */
  h_1 = 10223,
    /**< \f$\mbox{h}_{1}\f$
         (<code>h_1</code>) */
  Kstar_00 = 10311,
    /**< \f$\mbox{K}^{*0 }\f$
         (<code>K*_00</code>) */
  Kstar_0bar0 = -10311,
    /**< \f$\overline{\mbox{K}}^{*0 }\f$
         (<code>K*_0bar0</code>) */
  K_10 = 10313,
    /**< \f$\mbox{K}^{0 }_{1}\f$
         (<code>K_10</code>) */
  K_1bar0 = -10313,
    /**< \f$\overline{\mbox{K}}^{0 }_{1}\f$
         (<code>K_1bar0</code>) */
  Kstar_0plus = 10321,
    /**< \f$\mbox{K}^{*0 +}\f$
         (<code>K*_0+</code>) */
  Kstar_0minus = -10321,
    /**< \f$\mbox{K}^{*0 -}\f$
         (<code>K*_0-</code>) */
  K_1plus = 10323,
    /**< \f$\mbox{K}^{+}_{1}\f$
         (<code>K_1+</code>) */
  K_1minus = -10323,
    /**< \f$\mbox{K}^{-}_{1}\f$
         (<code>K_1-</code>) */
  eta1440 = 100331,
    /**< \f$\eta ^{0 }_{144}\f$
         (<code>eta1440</code>) */
  hprime_1 = 10333,
    /**< \f$\mbox{h}^{\prime }_{1}\f$
         (<code>h'_1</code>) */
  Dstar_0plus = 10411,
    /**< \f$\mbox{D}^{*0 +}\f$
         (<code>D*_0+</code>) */
  Dstar_0minus = -10411,
    /**< \f$\mbox{D}^{*0 -}\f$
         (<code>D*_0-</code>) */
  D_1plus = 10413,
    /**< \f$\mbox{D}^{+}_{1}\f$
         (<code>D_1+</code>) */
  D_1minus = -10413,
    /**< \f$\mbox{D}^{-}_{1}\f$
         (<code>D_1-</code>) */
  Dstar_00 = 10421,
    /**< \f$\mbox{D}^{*0 }\f$
         (<code>D*_00</code>) */
  Dstar_0bar0 = -10421,
    /**< \f$\overline{\mbox{D}}^{*0 }\f$
         (<code>D*_0bar0</code>) */
  D_10 = 10423,
    /**< \f$\mbox{D}^{0 }_{1}\f$
         (<code>D_10</code>) */
  D_1bar0 = -10423,
    /**< \f$\overline{\mbox{D}}^{0 }_{1}\f$
         (<code>D_1bar0</code>) */
  Dstar_0splus = 10431,
    /**< \f$\mbox{D}^{*+}_{0s}\f$
         (<code>D*_0s+</code>) */
  Dstar_0sminus = -10431,
    /**< \f$\mbox{D}^{*-}_{0s}\f$
         (<code>D*_0s-</code>) */
  D_1splus = 10433,
    /**< \f$\mbox{D}^{+}_{1s}\f$
         (<code>D_1s+</code>) */
  D_1sminus = -10433,
    /**< \f$\mbox{D}^{-}_{1s}\f$
         (<code>D_1s-</code>) */
  chi_0c = 10441,
    /**< \f$\chi _{0c}\f$
         (<code>chi_0c</code>) */
  h_1c = 10443,
    /**< \f$\mbox{h}_{1c}\f$
         (<code>h_1c</code>) */
  Bstar_00 = 10511,
    /**< \f$\mbox{B}^{*0 }\f$
         (<code>B*_00</code>) */
  Bstar_0bar0 = -10511,
    /**< \f$\overline{\mbox{B}}^{*0 }\f$
         (<code>B*_0bar0</code>) */
  B_10 = 10513,
    /**< \f$\mbox{B}^{0 }_{1}\f$
         (<code>B_10</code>) */
  B_1bar0 = -10513,
    /**< \f$\overline{\mbox{B}}^{0 }_{1}\f$
         (<code>B_1bar0</code>) */
  Bstar_0plus = 10521,
    /**< \f$\mbox{B}^{*0 +}\f$
         (<code>B*_0+</code>) */
  Bstar_0minus = -10521,
    /**< \f$\mbox{B}^{*0 -}\f$
         (<code>B*_0-</code>) */
  B_1plus = 10523,
    /**< \f$\mbox{B}^{+}_{1}\f$
         (<code>B_1+</code>) */
  B_1minus = -10523,
    /**< \f$\mbox{B}^{-}_{1}\f$
         (<code>B_1-</code>) */
  Bstar_0s0 = 10531,
    /**< \f$\mbox{B}^{*0 }_{0s}\f$
         (<code>B*_0s0</code>) */
  Bstar_0sbar0 = -10531,
    /**< \f$\overline{\mbox{B}}^{*0 }_{0s}\f$
         (<code>B*_0sbar0</code>) */
  B_1s0 = 10533,
    /**< \f$\mbox{B}^{0 }_{1s}\f$
         (<code>B_1s0</code>) */
  B_1sbar0 = -10533,
    /**< \f$\overline{\mbox{B}}^{0 }_{1s}\f$
         (<code>B_1sbar0</code>) */
  Bstar_0cplus = 10541,
    /**< \f$\mbox{B}^{*+}_{0c}\f$
         (<code>B*_0c+</code>) */
  Bstar_0cminus = -10541,
    /**< \f$\mbox{B}^{*-}_{0c}\f$
         (<code>B*_0c-</code>) */
  B_1cplus = 10543,
    /**< \f$\mbox{B}^{+}_{1c}\f$
         (<code>B_1c+</code>) */
  B_1cminus = -10543,
    /**< \f$\mbox{B}^{-}_{1c}\f$
         (<code>B_1c-</code>) */
  chi_0b = 10551,
    /**< \f$\chi _{0b}\f$
         (<code>chi_0b</code>) */
  h_1b = 10553,
    /**< \f$\mbox{h}_{1b}\f$
         (<code>h_1b</code>) */
  a_10 = 20113,
    /**< \f$\mbox{a}^{0 }_{1}\f$
         (<code>a_10</code>) */
  a_1plus = 20213,
    /**< \f$\mbox{a}^{+}_{1}\f$
         (<code>a_1+</code>) */
  a_1minus = -20213,
    /**< \f$\mbox{a}^{-}_{1}\f$
         (<code>a_1-</code>) */
  f_1 = 20223,
    /**< \f$\mbox{f}_{1}\f$
         (<code>f_1</code>) */
  Kstar_10 = 20313,
    /**< \f$\mbox{K}^{*0 }_{1}\f$
         (<code>K*_10</code>) */
  Kstar_1bar0 = -20313,
    /**< \f$\overline{\mbox{K}}^{*0 }_{1}\f$
         (<code>K*_1bar0</code>) */
  Kstar_1plus = 20323,
    /**< \f$\mbox{K}^{*+}_{1}\f$
         (<code>K*_1+</code>) */
  Kstar_1minus = -20323,
    /**< \f$\mbox{K}^{*-}_{1}\f$
         (<code>K*_1-</code>) */
  fprime_1 = 20333,
    /**< \f$\mbox{f}^{\prime }_{1}\f$
         (<code>f'_1</code>) */
  Dstar_1plus = 20413,
    /**< \f$\mbox{D}^{*+}_{1}\f$
         (<code>D*_1+</code>) */
  Dstar_1minus = -20413,
    /**< \f$\mbox{D}^{*-}_{1}\f$
         (<code>D*_1-</code>) */
  Dstar_10 = 20423,
    /**< \f$\mbox{D}^{*0 }_{1}\f$
         (<code>D*_10</code>) */
  Dstar_1bar0 = -20423,
    /**< \f$\overline{\mbox{D}}^{*0 }_{1}\f$
         (<code>D*_1bar0</code>) */
  Dstar_1splus = 20433,
    /**< \f$\mbox{D}^{*+}_{1s}\f$
         (<code>D*_1s+</code>) */
  Dstar_1sminus = -20433,
    /**< \f$\mbox{D}^{*-}_{1s}\f$
         (<code>D*_1s-</code>) */
  chi_1c = 20443,
    /**< \f$\chi _{1c}\f$
         (<code>chi_1c</code>) */
  Bstar_10 = 20513,
    /**< \f$\mbox{B}^{*0 }_{1}\f$
         (<code>B*_10</code>) */
  Bstar_1bar0 = -20513,
    /**< \f$\overline{\mbox{B}}^{*0 }_{1}\f$
         (<code>B*_1bar0</code>) */
  Bstar_1plus = 20523,
    /**< \f$\mbox{B}^{*+}_{1}\f$
         (<code>B*_1+</code>) */
  Bstar_1minus = -20523,
    /**< \f$\mbox{B}^{*-}_{1}\f$
         (<code>B*_1-</code>) */
  Bstar_1s0 = 20533,
    /**< \f$\mbox{B}^{*0 }_{1s}\f$
         (<code>B*_1s0</code>) */
  Bstar_1sbar0 = -20533,
    /**< \f$\overline{\mbox{B}}^{*0 }_{1s}\f$
         (<code>B*_1sbar0</code>) */
  Bstar_1cplus = 20543,
    /**< \f$\mbox{B}^{*+}_{1c}\f$
         (<code>B*_1c+</code>) */
  Bstar_1cminus = -20543,
    /**< \f$\mbox{B}^{*-}_{1c}\f$
         (<code>B*_1c-</code>) */
  chi_1b = 20553,
    /**< \f$\chi _{1b}\f$
         (<code>chi_1b</code>) */
  psiprime = 100443,
    /**< \f$\psi ^{\prime }\f$
         (<code>psi'</code>) */
  Upsilonprime = 100553,
    /**< \f$\Upsilon  ^{\prime }\f$
         (<code>Upsilon'</code>) */
  SUSY_d_L = 1000001,
    /**< \f$\tilde{\mbox{d}}_{L}\f$
         (<code>~d_L</code>) */
  SUSY_d_Lbar = -1000001,
    /**< \f$\tilde{\overline{\mbox{d}}}_{L}\f$
         (<code>~d_Lbar</code>) */
  SUSY_u_L = 1000002,
    /**< \f$\tilde{\mbox{u}}_{L}\f$
         (<code>~u_L</code>) */
  SUSY_u_Lbar = -1000002,
    /**< \f$\tilde{\overline{\mbox{u}}}_{L}\f$
         (<code>~u_Lbar</code>) */
  SUSY_s_L = 1000003,
    /**< \f$\tilde{\mbox{s}}_{L}\f$
         (<code>~s_L</code>) */
  SUSY_s_Lbar = -1000003,
    /**< \f$\tilde{\overline{\mbox{s}}}_{L}\f$
         (<code>~s_Lbar</code>) */
  SUSY_c_L = 1000004,
    /**< \f$\tilde{\mbox{c}}_{L}\f$
         (<code>~c_L</code>) */
  SUSY_c_Lbar = -1000004,
    /**< \f$\tilde{\overline{\mbox{c}}}_{L}\f$
         (<code>~c_Lbar</code>) */
  SUSY_b_1 = 1000005,
    /**< \f$\tilde{\mbox{b}}_{1}\f$
         (<code>~b_1</code>) */
  SUSY_b_1bar = -1000005,
    /**< \f$\tilde{\overline{\mbox{b}}}_{1}\f$
         (<code>~b_1bar</code>) */
  SUSY_t_1 = 1000006,
    /**< \f$\tilde{\mbox{t}}_{1}\f$
         (<code>~t_1</code>) */
  SUSY_t_1bar = -1000006,
    /**< \f$\tilde{\overline{\mbox{t}}}_{1}\f$
         (<code>~t_1bar</code>) */
  SUSY_e_Lminus = 1000011,
    /**< \f$\tilde{\mbox{e}}^{-}_{L}\f$
         (<code>~e_L-</code>) */
  SUSY_e_Lplus = -1000011,
    /**< \f$\tilde{\mbox{e}}^{+}_{L}\f$
         (<code>~e_L+</code>) */
  SUSY_nu_eL = 1000012,
    /**< \f$\tilde{\nu }_{eL}\f$
         (<code>~nu_eL</code>) */
  SUSY_nu_eLbar = -1000012,
    /**< \f$\tilde{\overline{\nu }}_{eL}\f$
         (<code>~nu_eLbar</code>) */
  SUSY_mu_Lminus = 1000013,
    /**< \f$\tilde{\mu }^{-}_{L}\f$
         (<code>~mu_L-</code>) */
  SUSY_mu_Lplus = -1000013,
    /**< \f$\tilde{\mu }^{+}_{L}\f$
         (<code>~mu_L+</code>) */
  SUSY_nu_muL = 1000014,
    /**< \f$\tilde{\nu }_{\mu L}\f$
         (<code>~nu_muL</code>) */
  SUSY_nu_muLbar = -1000014,
    /**< \f$\tilde{\overline{\nu }}_{\mu L}\f$
         (<code>~nu_muLbar</code>) */
  SUSY_tau_1minus = 1000015,
    /**< \f$\tilde{\tau }^{-}_{1}\f$
         (<code>~tau_1-</code>) */
  SUSY_tau_1plus = -1000015,
    /**< \f$\tilde{\tau }^{+}_{1}\f$
         (<code>~tau_1+</code>) */
  SUSY_nu_tauL = 1000016,
    /**< \f$\tilde{\nu }_{\tau L}\f$
         (<code>~nu_tauL</code>) */
  SUSY_nu_tauLbar = -1000016,
    /**< \f$\tilde{\overline{\nu }}_{\tau L}\f$
         (<code>~nu_tauLbar</code>) */
  SUSY_g = 1000021,
    /**< \f$\tilde{\mbox{g}}\f$
         (<code>~g</code>) */
  SUSY_chi_10 = 1000022,
    /**< \f$\tilde{\chi }^{0 }_{1}\f$
         (<code>~chi_10</code>) */
  SUSY_chi_20 = 1000023,
    /**< \f$\tilde{\chi }^{0 }_{2}\f$
         (<code>~chi_20</code>) */
  SUSY_chi_1plus = 1000024,
    /**< \f$\tilde{\chi }^{+}_{1}\f$
         (<code>~chi_1+</code>) */
  SUSY_chi_1minus = -1000024,
    /**< \f$\tilde{\chi }^{-}_{1}\f$
         (<code>~chi_1-</code>) */
  SUSY_chi_30 = 1000025,
    /**< \f$\tilde{\chi }^{0 }_{3}\f$
         (<code>~chi_30</code>) */
  SUSY_chi_40 = 1000035,
    /**< \f$\tilde{\chi }^{0 }_{4}\f$
         (<code>~chi_40</code>) */
  SUSY_chi_2plus = 1000037,
    /**< \f$\tilde{\chi }^{+}_{2}\f$
         (<code>~chi_2+</code>) */
  SUSY_chi_2minus = -1000037,
    /**< \f$\tilde{\chi }^{-}_{2}\f$
         (<code>~chi_2-</code>) */
  SUSY_Gravitino = 1000039,
    /**< \f$\tilde{{\cal G}}\f$
         (<code>~Gravitino</code>) */
  SUSY_d_R = 2000001,
    /**< \f$\tilde{\mbox{d}}_{R}\f$
         (<code>~d_R</code>) */
  SUSY_d_Rbar = -2000001,
    /**< \f$\tilde{\overline{\mbox{d}}}_{R}\f$
         (<code>~d_Rbar</code>) */
  SUSY_u_R = 2000002,
    /**< \f$\tilde{\mbox{u}}_{R}\f$
         (<code>~u_R</code>) */
  SUSY_u_Rbar = -2000002,
    /**< \f$\tilde{\overline{\mbox{u}}}_{R}\f$
         (<code>~u_Rbar</code>) */
  SUSY_s_R = 2000003,
    /**< \f$\tilde{\mbox{s}}_{R}\f$
         (<code>~s_R</code>) */
  SUSY_s_Rbar = -2000003,
    /**< \f$\tilde{\overline{\mbox{s}}}_{R}\f$
         (<code>~s_Rbar</code>) */
  SUSY_c_R = 2000004,
    /**< \f$\tilde{\mbox{c}}_{R}\f$
         (<code>~c_R</code>) */
  SUSY_c_Rbar = -2000004,
    /**< \f$\tilde{\overline{\mbox{c}}}_{R}\f$
         (<code>~c_Rbar</code>) */
  SUSY_b_2 = 2000005,
    /**< \f$\tilde{\mbox{b}}_{2}\f$
         (<code>~b_2</code>) */
  SUSY_b_2bar = -2000005,
    /**< \f$\tilde{\overline{\mbox{b}}}_{2}\f$
         (<code>~b_2bar</code>) */
  SUSY_t_2 = 2000006,
    /**< \f$\tilde{\mbox{t}}_{2}\f$
         (<code>~t_2</code>) */
  SUSY_t_2bar = -2000006,
    /**< \f$\tilde{\overline{\mbox{t}}}_{2}\f$
         (<code>~t_2bar</code>) */
  SUSY_e_Rminus = 2000011,
    /**< \f$\tilde{\mbox{e}}^{-}_{R}\f$
         (<code>~e_R-</code>) */
  SUSY_e_Rplus = -2000011,
    /**< \f$\tilde{\mbox{e}}^{+}_{R}\f$
         (<code>~e_R+</code>) */
  SUSY_nu_eR = 2000012,
    /**< \f$\tilde{\nu }_{eR}\f$
         (<code>~nu_eR</code>) */
  SUSY_nu_eRbar = -2000012,
    /**< \f$\tilde{\overline{\nu }}_{eR}\f$
         (<code>~nu_eRbar</code>) */
  SUSY_mu_Rminus = 2000013,
    /**< \f$\tilde{\mu }^{-}_{R}\f$
         (<code>~mu_R-</code>) */
  SUSY_mu_Rplus = -2000013,
    /**< \f$\tilde{\mu }^{+}_{R}\f$
         (<code>~mu_R+</code>) */
  SUSY_nu_muR = 2000014,
    /**< \f$\tilde{\nu }_{\mu R}\f$
         (<code>~nu_muR</code>) */
  SUSY_nu_muRbar = -2000014,
    /**< \f$\tilde{\overline{\nu }}_{\mu R}\f$
         (<code>~nu_muRbar</code>) */
  SUSY_tau_2minus = 2000015,
    /**< \f$\tilde{\tau }^{-}_{2}\f$
         (<code>~tau_2-</code>) */
  SUSY_tau_2plus = -2000015,
    /**< \f$\tilde{\tau }^{+}_{2}\f$
         (<code>~tau_2+</code>) */
  SUSY_nu_tauR = 2000016,
    /**< \f$\tilde{\nu }_{\tau R}\f$
         (<code>~nu_tauR</code>) */
  SUSY_nu_tauRbar = -2000016,
    /**< \f$\tilde{\overline{\nu }}_{\tau R}\f$
         (<code>~nu_tauRbar</code>) */
  pi_tc0 = 3000111,
    /**< \f$\pi ^{0 }_{tc}\f$
         (<code>pi_tc0</code>) */
  pi_tcplus = 3000211,
    /**< \f$\pi ^{+}_{tc}\f$
         (<code>pi_tc+</code>) */
  pi_tcminus = -3000211,
    /**< \f$\pi ^{-}_{tc}\f$
         (<code>pi_tc-</code>) */
  piprime_tc0 = 3000221,
    /**< \f$\pi ^{\prime 0 }_{tc}\f$
         (<code>pi'_tc0</code>) */
  eta_tc0 = 3000331,
    /**< \f$\eta ^{0 }_{tc}\f$
         (<code>eta_tc0</code>) */
  rho_tc0 = 3000113,
    /**< \f$\rho ^{0 }_{tc}\f$
         (<code>rho_tc0</code>) */
  rho_tcplus = 3000213,
    /**< \f$\rho ^{+}_{tc}\f$
         (<code>rho_tc+</code>) */
  rho_tcminus = -3000213,
    /**< \f$\rho ^{-}_{tc}\f$
         (<code>rho_tc-</code>) */
  omega_tc = 3000223,
    /**< \f$\omega _{tc}\f$
         (<code>omega_tc</code>) */
  V8_tc = 3100021,
    /**< \f$\mbox{V}_{8tc}\f$
         (<code>V8_tc</code>) */
  pi_22_1_tc = 3100111,
    /**< \f$\pi _{22}\f$
         (<code>pi_22_1_tc</code>) */
  pi_22_8_tc = 3200111,
    /**< \f$\pi _{22}\f$
         (<code>pi_22_8_tc</code>) */
  rho_11_tc = 3100113,
    /**< \f$\rho _{11}\f$
         (<code>rho_11_tc</code>) */
  rho_12_tc = 3200113,
    /**< \f$\rho _{12}\f$
         (<code>rho_12_tc</code>) */
  rho_21_tc = 3300113,
    /**< \f$\rho _{21}\f$
         (<code>rho_21_tc</code>) */
  rho_22_tc = 3400113,
    /**< \f$\rho _{22}\f$
         (<code>rho_22_tc</code>) */
  dstar = 4000001,
    /**< \f$\mbox{d}^{*}\f$
         (<code>d*</code>) */
  dstarbar = -4000001,
    /**< \f$\overline{\mbox{d}}^{*}\f$
         (<code>d*bar</code>) */
  ustar = 4000002,
    /**< \f$\mbox{u}^{*}\f$
         (<code>u*</code>) */
  ustarbar = -4000002,
    /**< \f$\overline{\mbox{u}}^{*}\f$
         (<code>u*bar</code>) */
  estarminus = 4000011,
    /**< \f$\mbox{e}^{*-}\f$
         (<code>e*-</code>) */
  estarbarplus = -4000011,
    /**< \f$\overline{\mbox{e}}^{*+}\f$
         (<code>e*bar+</code>) */
  nustar_e0 = 4000012,
    /**< \f$\nu ^{*0 }_{e}\f$
         (<code>nu*_e0</code>) */
  nustar_ebar0 = -4000012,
    /**< \f$\overline{\nu }^{*0 }_{e}\f$
         (<code>nu*_ebar0</code>) */
  Gravitonstar = 5000039,
    /**< \f${\cal G}^{*}\f$
         (<code>Graviton*</code>) */
  nu_Re = 9900012,
    /**< \f$\nu _{Re}\f$
         (<code>nu_Re</code>) */
  nu_Rmu = 9900014,
    /**< \f$\nu _{R\mu }\f$
         (<code>nu_Rmu</code>) */
  nu_Rtau = 9900016,
    /**< \f$\nu _{R\tau }\f$
         (<code>nu_Rtau</code>) */
  Z_R0 = 9900023,
    /**< \f$\mbox{Z}^{0 }_{R}\f$
         (<code>Z_R0</code>) */
  W_Rplus = 9900024,
    /**< \f$\mbox{W}^{+}_{R}\f$
         (<code>W_R+</code>) */
  W_Rminus = -9900024,
    /**< \f$\mbox{W}^{-}_{R}\f$
         (<code>W_R-</code>) */
  H_Lplus2 = 9900041,
    /**< \f$\mbox{H}^{++}_{L}\f$
         (<code>H_L++</code>) */
  H_Lminus2 = -9900041,
    /**< \f$\mbox{H}^{--}_{L}\f$
         (<code>H_L--</code>) */
  H_Rplus2 = 9900042,
    /**< \f$\mbox{H}^{++}_{R}\f$
         (<code>H_R++</code>) */
  H_Rminus2 = -9900042,
    /**< \f$\mbox{H}^{--}_{R}\f$
         (<code>H_R--</code>) */
  rho_diff0 = 9900110,
    /**< \f$\rho ^{0 }_{\mbox{\scriptsize diffr}}\f$
         (<code>rho_diff0</code>) */
  pi_diffrplus = 9900210,
    /**< \f$\pi ^{+}_{\mbox{\scriptsize diffr}}\f$
         (<code>pi_diffr+</code>) */
  pi_diffrminus = -9900210,
    /**< \f$\pi ^{-}_{\mbox{\scriptsize diffr}}\f$
         (<code>pi_diffr-</code>) */
  omega_di = 9900220,
    /**< \f$\omega _{\mbox{\scriptsize diffr}}\f$
         (<code>omega_di</code>) */
  phi_diff = 9900330,
    /**< \f$\phi _{\mbox{\scriptsize diffr}}\f$
         (<code>phi_diff</code>) */
  Jpsi_di = 9900440,
    /**< \f$J/\psi _{\mbox{\scriptsize diffr}}\f$
         (<code>J/psi_di</code>) */
  n_diffr0 = 9902110,
    /**< \f$\mbox{n}^{0 }_{\mbox{\scriptsize diffr}}\f$
         (<code>n_diffr0</code>) */
  n_diffrbar0 = -9902110,
    /**< \f$\overline{\mbox{n}}^{0 }_{\mbox{\scriptsize diffr}}\f$
         (<code>n_diffrbar0</code>) */
  p_diffrplus = 9902210,
    /**< \f$\mbox{p}^{+}_{\mbox{\scriptsize diffr}}\f$
         (<code>p_diffr+</code>) */
  p_diffrbarminus = -9902210,
    /**< \f$\overline{\mbox{p}}^{-}_{\mbox{\scriptsize diffr}}\f$
         (<code>p_diffrbar-</code>) */
  undefined = 0
    /**< Undefined particle. */ };
}
}
#endif
