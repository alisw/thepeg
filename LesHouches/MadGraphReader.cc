// -*- C++ -*-
//
// MadGraphReader.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MadGraphReader class.
//

#include "MadGraphReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/PDFBase.h"
#include "ThePEG/LesHouches/MadGraphOneCut.h"
#include "ThePEG/LesHouches/MadGraphTwoCut.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using std::fgetc;
using std::fgets;

IBPtr MadGraphReader::clone() const {
  return new_ptr(*this);
}

IBPtr MadGraphReader::fullclone() const {
  return new_ptr(*this);
}

void MadGraphReader::open() {
  LesHouchesFileReader::open();

  heprup.IDWTUP = 1;

  static const int ntags = 33;
  static const char * cuttags[] = {"ptj", "ptb", "pta", "ptl",
				   "etaj", "etab", "etaa", "etal",
				   "drjj", "drbb", "draa", "drll",
				   "drbj", "draj", "drjl", "drab",
				   "drbl", "dral", "mmjj", "mmbb",
				   "mmaa", "mmll", "mmbj", "mmaj",
				   "mmjl", "mmab", "mmbl", "mmal", 
				   "xptj", "xptb", "xpta", "xptl", "xetamin"};

  ieve = neve = 0;

  // If we are reading a LHF formatted file things are rather easy.

  if ( LHFVersion.size() ) {
    // extract the number of events
    neve = numberOfEvents(outsideBlock);
    string madheader = outsideBlock;
    if ( neve == 0 ) {
      neve = numberOfEvents(headerBlock);
      madheader = headerBlock;
    }
    if ( neve == 0 )
      Throw<LesHouchesFileError>()
      << "The MadGraphReader '" << name() << "' expected the LHE file '"
      << filename() << "' to include MadGraph-specific header information, "
      << "but did not find any. The events may not be properly sampled."
      << Exception::warning;

    if ( neve != 0 ) NEvents(neve);

    // MadEvent has gives wrong values for XMAXUP and XWGTUP, they
    // need to be multiplied by the number of events to be LHEF
    // compliant.
    weightScale = neve*picobarn;

    // Scan information about cuts.
    for ( int itag = 0; itag < ntags; ++itag ) {
     string::size_type pos = madheader.find(string("= ") + cuttags[itag]);
      if ( pos != string::npos ) {
	string::size_type beg = max(madheader.rfind("#", pos) + 1,
				    madheader.rfind("\n", pos) + 1);
	string value = madheader.substr(beg, pos - beg);
	for ( string::size_type i = 0; i < value.length(); ++i )
	  if ( value[i] == 'd' || value[i] == 'D' ) value[i] = 'e';
	cuts[cuttags[itag]] = std::strtod(value.c_str(), NULL);
      }
    }

    return;

  }

  double xsec = -1.0;
  double maxw = -1.0;
  double ebeam1 = -1.0;
  double ebeam2 = -1.0;
  int lpp1 = 0;
  int lpp2 = 0;
  string pdftag;
  // First scan banner to extract some information
  // (LesHoushesFileReader::open has already read in the first line).
  cfile.resetline();
  do {
    if ( !cfile ) break;
    if ( cfile.getc() != '#' ) {
      int test;
      cfile >> test;
      if ( cfile ) {
	cfile.resetline();
	break;
      }
    }
    if ( cfile.find("#  Number of Events") ) {
      cfile.skip(':');
      cfile >> neve;
    } else if ( cfile.find("Integrated weight") ) {
      cfile.skip(':');
      cfile >> xsec;
    } else if ( cfile.find("Max wgt") ) {
      cfile.skip(':');
      cfile >> maxw;
    } else if ( cfile.find("ebeam(1)") || cfile.find("ebeam1") ) {
      cfile >> ebeam1;
    } else if ( cfile.find("ebeam(2)") || cfile.find("ebeam2") ) {
      cfile >> ebeam2;
    } else if ( cfile.find("lpp(1)") || cfile.find("lpp1") ) {
      cfile >> lpp1;
    } else if ( cfile.find("lpp(2)") || cfile.find("lpp2") ) {
      cfile >> lpp2;
    } else if ( cfile.find("PDF set") ) {
      cfile.skip('\'');
      cfile >> pdftag;
      pdftag = pdftag.substr(0, 7);
    } else if ( cfile.find("Number of Events Written ") ) {
      cfile.skip(':');
      cfile >> neve;
      maxw = xsec/double(neve);
    } else {
      for ( int itag = 0; itag < ntags; ++itag ) {
	if ( cfile.find(string("= ") + cuttags[itag]) ) {
	  
	  cfile >> cuts[cuttags[itag]];
	  if ( cfile.getc() == 'd' ) {
	    long x = 0;
	    cfile >> x;
	    cuts[cuttags[itag]] *= pow(10.0, double(x));
	  }
	  break;
	}
      }
    }
  } while ( cfile.readline() );

  // Return here if no comment block was found.
  if ( neve <= 0 ) return;

  // Convert the extracted information to LesHouches format.
  heprup.NPRUP = 1;
  heprup.LPRUP.push_back(0);
  heprup.XSECUP.push_back(xsec);
  heprup.XERRUP.push_back(0.0);
  heprup.XMAXUP.push_back(maxw);
  NEvents(neve);
  // MadEvent has gives wrong values for XMAXUP and XWGTUP, they
  // need to be multiplied by the number of events to be LHEF
  // compliant.
  weightScale = neve*picobarn;

  if ( !heprup.IDBMUP.first ) {
    if ( lpp1 == 1 ) heprup.IDBMUP.first = ParticleID::pplus;
    else if ( lpp1 == -1 ) heprup.IDBMUP.first = ParticleID::pbarminus;
  }
  if ( !heprup.IDBMUP.second ) {
    if ( lpp2 == 1 ) heprup.IDBMUP.second = ParticleID::pplus;
    else if ( lpp2 == -1 ) heprup.IDBMUP.second = ParticleID::pbarminus;
  }

  if ( heprup.EBMUP.first <= 0.0 ) heprup.EBMUP.first = ebeam1;
  if ( heprup.EBMUP.second <= 0.0 ) heprup.EBMUP.second = ebeam2;

  if ( !cfile )
    throw LesHouchesFileError()
      << "An error occurred while '" << name() << "' was reading the file '"
      << filename() << "'." << Exception::runerror;

  if ( heprup.PDFSUP.first != 0 || heprup.PDFSUP.first != 0 ) return;
  // If we have an old MadGraph we have to try to figure out which PDF
  // codes to use.
  heprup.PDFGUP.first = heprup.PDFGUP.second = 0;
  if ( pdftag == "mrs02nl" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 20200;
  else if ( pdftag == "mrs02nn" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 20270;
  else if ( pdftag == "cteq6_m" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 10050;
  else if ( pdftag == "cteq6_l" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 10041;
  else if ( pdftag == "cteq6l1" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 10042;
  else if ( pdftag == "cteq5_m" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 19050;
  else if ( pdftag == "cteq5_d" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 19060;
  else if ( pdftag == "cteq5_l" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 19070;
  else if ( pdftag == "cteq4_m" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 19150;
  else if ( pdftag == "cteq4_d" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 19160;
  else if ( pdftag == "cteq4_l" )
    heprup.PDFSUP.first = heprup.PDFSUP.second = 19170;
}


long MadGraphReader::scan() {
  bool fixscale = !NEvents();
  long neve = LesHouchesFileReader::scan();
  if ( fixscale ) {
    // MadEvent has gives wrong values for XMAXUP and XWGTUP, they
    // need to be multiplied by the number of events to be LHEF
    // compliant.
    weightScale = neve*picobarn;
    if ( heprup.NPRUP > 1 ) weightScale /= heprup.NPRUP;
  }
  return neve;
}

bool MadGraphReader::doReadEvent() {

  if ( LesHouchesFileReader::doReadEvent() )  return true;

  if ( !cfile ) return false;

  hepeup.NUP = 0;
  ieve = 0;
  long evno = 0;
  hepeup.XWGTUP = 0.0;
  double scale = 0.0;
  double aEM = 0.0;
  double aS = 0.0;
  bool oldformat = false;

  cfile >> hepeup.NUP >> evno >> hepeup.XWGTUP >> scale >> aEM >> aS;
  if ( !cfile ) {
    hepeup.IDPRUP = evno;
    hepeup.SCALUP = fixedScale/GeV;
    hepeup.AQEDUP = fixedAEM;
    hepeup.AQCDUP = fixedAS;
    ++ieve;
    oldformat = true;
  } else {
    hepeup.IDPRUP = 0;
    ieve = evno;
    hepeup.SCALUP = scale;
    hepeup.AQEDUP = aEM;
    hepeup.AQCDUP = aS;
  }

  hepeup.IDUP.resize(hepeup.NUP);
  if ( !cfile.readline() ) return false;
  for ( int i = 0; i < hepeup.NUP; ++i ) cfile >> hepeup.IDUP[i];

  hepeup.MOTHUP.resize(hepeup.NUP);
  if ( !cfile.readline() ) return false;
  for ( int i = 0; i < hepeup.NUP; ++i ) cfile >> hepeup.MOTHUP[i].first;
  if ( !cfile.readline() ) return false;
  for ( int i = 0; i < hepeup.NUP; ++i ) cfile >> hepeup.MOTHUP[i].second;

  hepeup.ICOLUP.resize(hepeup.NUP);
  if ( !cfile.readline() ) return false;
  for ( int i = 0; i < hepeup.NUP; ++i ) cfile >> hepeup.ICOLUP[i].first;
  if ( !cfile.readline() ) return false;
  for ( int i = 0; i < hepeup.NUP; ++i ) cfile >> hepeup.ICOLUP[i].second;

  // Try to figure out if the colour lines are reversed
  bool colrev = false;
  for ( int i = 0; i < hepeup.NUP; ++i )
    if ( abs(hepeup.IDUP[i]) < 10 && hepeup.IDUP[i] < 0 &&
	 !hepeup.ICOLUP[i].second ) colrev = true;
  if ( colrev ) for ( int i = 0; i < hepeup.NUP; ++i )
    swap(hepeup.ICOLUP[i].first, hepeup.ICOLUP[i].second);

  if ( oldformat ) {
    hepeup.ISTUP.assign(hepeup.NUP, 1);
    hepeup.ISTUP[0] = hepeup.ISTUP[1] = -1;
    hepeup.SPINUP.assign(hepeup.NUP, 9);
  } else {
    hepeup.ISTUP.resize(hepeup.NUP);
    if ( !cfile.readline() ) return false;
    for ( int i = 0; i < hepeup.NUP; ++i ) cfile >> hepeup.ISTUP[i];
    hepeup.SPINUP.resize(hepeup.NUP, 9);
    if ( !cfile.readline() ) return false;
    for ( int i = 0; i < hepeup.NUP; ++i ) cfile >> hepeup.SPINUP[i];
  }

  hepeup.PUP.resize(hepeup.NUP, vector<double>(5));
  for ( int i = 0; i < hepeup.NUP; ++i ) {
    if ( !cfile.readline() ) return false;
    int dummy = 0;
    cfile  >> dummy >> hepeup.PUP[i][3]
	   >> hepeup.PUP[i][0] >> hepeup.PUP[i][1] >> hepeup.PUP[i][2];
    hepeup.PUP[i][4] =
      sqrt(max(sqr(hepeup.PUP[i][3]) - sqr(hepeup.PUP[i][0]) -
	       sqr(hepeup.PUP[i][1]) - sqr(hepeup.PUP[i][2]), 0.0));
  }

  if ( !cfile ) return false;

  // Set info not obtained from MadGraph.
  hepeup.VTIMUP = vector<double>(hepeup.NUP, -1.0);

  // Deduce positions of incoming beams and corresponding partons.
  pair<int,int> beampos(-1, -1);
  for ( int i = 0; i < hepeup.NUP; ++i ) {
    if ( hepeup.ISTUP[i] != -9 ) continue;
    if ( beampos.first < 0 ) beampos.first = i;
    else if ( beampos.second < 0 ) beampos.second = i;
  }
  pair<int,int> partpos(-1, -1);
  for ( int i = hepeup.NUP - 1; i >= 0; --i ) {
    if ( hepeup.ISTUP[i] != -1 ) continue;
    if ( hepeup.MOTHUP[i].first > 0 &&
	 hepeup.MOTHUP[i].first == beampos.first )
      partpos.first = i;
    else if ( hepeup.MOTHUP[i].first > 0 &&
	      hepeup.MOTHUP[i].first == beampos.second )
      partpos.second = i;
    else if ( partpos.second < 0 ) partpos.second = i;
    else if ( partpos.first < 0 ) partpos.first = i;
  }

  // We set these to -1 to let the base class do the work.
  hepeup.XPDWUP.first = -1.0;
  hepeup.XPDWUP.second = -1.0;

  cfile.readline();

  // Return true even if last read failed.
  return true;

}

void MadGraphReader::persistentOutput(PersistentOStream & os) const {
  os << ounit(fixedScale, GeV) << fixedAEM << fixedAS << cuts
     << doInitCuts;
}

void MadGraphReader::persistentInput(PersistentIStream & is, int) {
  is >> iunit(fixedScale, GeV) >> fixedAEM >> fixedAS >> cuts
     >> doInitCuts;
}

bool MadGraphReader::preInitialize() const {
  if ( LesHouchesFileReader::preInitialize() ) return true;
  if ( doInitCuts && !theCuts ) return true;
  return false;
}

void MadGraphReader::doinit() {
  LesHouchesFileReader::doinit();
  if ( doInitCuts && !theCuts ) {
    theCuts = initCuts();
    if ( !theCuts ) Throw<Exception>()
      << "MadGraphReader '" << name()
      << "' could not create cut objects in pre-initialization."
      << Exception::warning;
  }
}

void MadGraphReader::initPDFs() {
  LesHouchesFileReader::initPDFs();
}

CutsPtr MadGraphReader::initCuts() {
  CutsPtr newCuts;
  open();
  close();
  if ( cuts.empty() ) return CutsPtr();
  vector<OneCutPtr> ones;
  vector<TwoCutPtr> twos;
  vector<string> onames;
  vector<string> tnames;
  for ( map<string,double>::iterator i = cuts.begin(); i != cuts.end(); ++i ) {
    if ( i->second <= 0.0 ) continue;
    MadGraphOneCut::CutType t = MadGraphOneCut::PT;
    char p = 0;
    if ( i->first.substr(0, 2) == "pt" ) {
      t = MadGraphOneCut::PT;
      p = i->first[2];
    }
    else if  ( i->first.substr(0, 3) == "eta" ) {
      t = MadGraphOneCut::ETA;
      p = i->first[3];
    }
    else if  ( i->first.substr(0, 3) == "xpt" ) {
      t = MadGraphOneCut::XPT;
      p = i->first[3];
    }
    if ( p ) {
      MadGraphOneCut::PType pt = MadGraphOneCut::JET;
      switch ( p ) {
      case 'j':	pt = MadGraphOneCut::JET; break;
      case 'b':	pt = MadGraphOneCut::BOT; break;
      case 'a':	pt = MadGraphOneCut::PHO; break;
      case 'l':	pt = MadGraphOneCut::LEP; break;
      }
      ones.push_back(new_ptr(MadGraphOneCut(t, pt, i->second)));
      onames.push_back(i->first);
      continue;
    }
    if ( i->first.substr(0, 2) == "dr" || i->first.substr(0, 2) == "mm" ) {
      MadGraphTwoCut::CutType tt = MadGraphTwoCut::DELTAR;
      if ( i->first.substr(0, 2) == "mm" ) tt = MadGraphTwoCut::INVMASS;
      MadGraphTwoCut::PPType pp = MadGraphTwoCut::JETJET;
      if ( i->first.substr(2, 2) == "jj" ) pp = MadGraphTwoCut::JETJET;
      else if ( i->first.substr(2, 2) == "bb" )	pp = MadGraphTwoCut::BOTBOT;
      else if ( i->first.substr(2, 2) == "aa" )	pp = MadGraphTwoCut::PHOPHO;
      else if ( i->first.substr(2, 2) == "ll" )	pp = MadGraphTwoCut::LEPLEP;
      else if ( i->first.substr(2, 2) == "bj" )	pp = MadGraphTwoCut::BOTJET;
      else if ( i->first.substr(2, 2) == "aj" )	pp = MadGraphTwoCut::PHOJET;
      else if ( i->first.substr(2, 2) == "jl" )	pp = MadGraphTwoCut::JETLEP;
      else if ( i->first.substr(2, 2) == "ab" )	pp = MadGraphTwoCut::PHOBOT;
      else if ( i->first.substr(2, 2) == "bl" )	pp = MadGraphTwoCut::BOTLEP;
      else if ( i->first.substr(2, 2) == "al" )	pp = MadGraphTwoCut::PHOLEP;
      twos.push_back(new_ptr(MadGraphTwoCut(tt, pp, i->second)));
      tnames.push_back(i->first);
    }
  }
  if ( ones.empty() && twos.empty() ) return CutsPtr();

  newCuts = new_ptr(Cuts());
  generator()->preinitRegister(newCuts, fullName() + "/ExtractedCuts");

  for ( int i = 0, N = ones.size(); i < N; ++i ) {
    generator()->preinitRegister(ones[i], fullName() + "/" + onames[i]);
    generator()->preinitInterface
      (newCuts, "OneCuts", 0, "insert",  ones[i]->fullName());
    //    newCuts->add(tOneCutPtr(ones[i]));
  }

  for ( int i = 0, N = twos.size(); i < N; ++i ) {
    reporeg(twos[i], tnames[i]);
    generator()->preinitInterface
      (newCuts, "TwoCuts", 0, "insert",  twos[i]->fullName());
    //    newCuts->add(tTwoCutPtr(twos[i]));
  }

  return newCuts;


}

string MadGraphReader::scanCuts(string) {
  if ( theCuts )
    return "A Cuts object has already been assigned to this reader.";
  open();
  close();
  if ( cuts.empty() ) return "No information about cuts were found. "
			"Maybe the file was from an old version of MadGraph";
  vector<OneCutPtr> ones;
  vector<TwoCutPtr> twos;
  vector<string> onames;
  vector<string> tnames;
  for ( map<string,double>::iterator i = cuts.begin(); i != cuts.end(); ++i ) {
    if ( i->second <= 0.0 ) continue;
    MadGraphOneCut::CutType t = MadGraphOneCut::PT;
    char p = 0;
    if ( i->first.substr(0, 2) == "pt" ) {
      t = MadGraphOneCut::PT;
      p = i->first[2];
    }
    else if  ( i->first.substr(0, 3) == "eta" ) {
      t = MadGraphOneCut::ETA;
      p = i->first[3];
    }
    else if  ( i->first.substr(0, 3) == "xpt" ) {
      t = MadGraphOneCut::XPT;
      p = i->first[3];
    }
    if ( p ) {
      MadGraphOneCut::PType pt = MadGraphOneCut::JET;
      switch ( p ) {
      case 'j':	pt = MadGraphOneCut::JET; break;
      case 'b':	pt = MadGraphOneCut::BOT; break;
      case 'a':	pt = MadGraphOneCut::PHO; break;
      case 'l':	pt = MadGraphOneCut::LEP; break;
      }
      ones.push_back(new_ptr(MadGraphOneCut(t, pt, i->second)));
      onames.push_back(i->first);
      continue;
    }
    if ( i->first.substr(0, 2) == "dr" || i->first.substr(0, 2) == "mm" ) {
      MadGraphTwoCut::CutType tt = MadGraphTwoCut::DELTAR;
      if ( i->first.substr(0, 2) == "mm" ) tt = MadGraphTwoCut::INVMASS;
      MadGraphTwoCut::PPType pp = MadGraphTwoCut::JETJET;
      if ( i->first.substr(2, 2) == "jj" ) pp = MadGraphTwoCut::JETJET;
      else if ( i->first.substr(2, 2) == "bb" )	pp = MadGraphTwoCut::BOTBOT;
      else if ( i->first.substr(2, 2) == "aa" )	pp = MadGraphTwoCut::PHOPHO;
      else if ( i->first.substr(2, 2) == "ll" )	pp = MadGraphTwoCut::LEPLEP;
      else if ( i->first.substr(2, 2) == "bj" )	pp = MadGraphTwoCut::BOTJET;
      else if ( i->first.substr(2, 2) == "aj" )	pp = MadGraphTwoCut::PHOJET;
      else if ( i->first.substr(2, 2) == "jl" )	pp = MadGraphTwoCut::JETLEP;
      else if ( i->first.substr(2, 2) == "ab" )	pp = MadGraphTwoCut::PHOBOT;
      else if ( i->first.substr(2, 2) == "bl" )	pp = MadGraphTwoCut::BOTLEP;
      else if ( i->first.substr(2, 2) == "al" )	pp = MadGraphTwoCut::PHOLEP;
      twos.push_back(new_ptr(MadGraphTwoCut(tt, pp, i->second)));
      tnames.push_back(i->first);
    }
  }
  if ( ones.empty() && twos.empty() ) return "No non-zero cuts found.";

  theCuts = new_ptr(Cuts());
  reporeg(theCuts, "ExtractedCuts");

  for ( int i = 0, N = ones.size(); i < N; ++i ) {
    reporeg(ones[i], onames[i]);
    theCuts->add(tOneCutPtr(ones[i]));
  }

  for ( int i = 0, N = twos.size(); i < N; ++i ) {
    reporeg(twos[i], tnames[i]);
    theCuts->add(tTwoCutPtr(twos[i]));
  }

  return "";
}

ClassDescription<MadGraphReader> MadGraphReader::initMadGraphReader;
// Definition of the static class description member.

void MadGraphReader::Init() {

  static ClassDocumentation<MadGraphReader> documentation
    ("ThePEG::MadGraphReader is used together with the LesHouchesEventHandler "
     "to read event files generated with the MadGraph/MadEvent program.",
     "Events were read from event files generated "
     "with the MadGraph/MadEvent\\cite{ThePEG::MadGraph} program.",
     "\\bibitem{ThePEG::MadGraph} F. Maltoni and T. Stelzer, "
     "hep-ph/0208156;\\\\"
     "T. Stelzer and W.F. Long, \\textit{Comput.~Phys.~Commun.} "
     "\\textbf{81} (1994) 357-371.");

  static Parameter<MadGraphReader,Energy> interfaceFixedScale
    ("FixedScale",
     "Old MadGraph files do not necessarily contain information about "
     "the factorization (or renormalization) scale. In this case this "
     "is used instead.",
     &MadGraphReader::fixedScale, GeV, 91.188*GeV, ZERO, 1000.0*GeV,
     true, false, true);
  interfaceFixedScale.setHasDefault(false);

  static Parameter<MadGraphReader,double> interfaceFixedAlphaEM
    ("FixedAlphaEM",
     "Old MadGraph files do not necessarily contain information about "
     "the value of \\f$\\alpha_{EM}\\f$. In this case this is used instead.",
     &MadGraphReader::fixedAEM, 0.007546772, 0.0, 1.0,
     true, false, true);
  interfaceFixedAlphaEM.setHasDefault(false);

  static Parameter<MadGraphReader,double> interfaceFixedAlphaS
    ("FixedAlphaS",
     "Old MadGraph files do not necessarily contain information about "
     "the value of \\f$\\alpha_S\\f$. In this case this is used instead.",
     &MadGraphReader::fixedAS, 0.12, 0.0, 1.0,
     true, false, true);
  interfaceFixedAlphaS.setHasDefault(false);

  static Command<MadGraphReader> interfaceScanCuts
    ("ScanCuts",
     "If no <interface>LesHouchesReader::Cuts</interface> has been assigned, "
     "the event file is scanned for information about generation cuts. If cuts "
     "are found, the corresponding objects will be created in a sub-directory "
     "with the same name as this object and assigned as the "
     "<interface>LesHouchesReader::Cuts</interface> of this reader.",
     &MadGraphReader::scanCuts, true);

  static Switch<MadGraphReader,bool> interfaceInitCuts
    ("InitCuts",
     "If no cuts were specified for this reader, try to extract cut "
     "information from the MadGraph file and assign the relevant cut "
     "objects when the reader is initialized.",
     &MadGraphReader::doInitCuts, false, true, false);
  static SwitchOption interfaceInitCutsYes
    (interfaceInitCuts,
     "Yes",
     "Extract cuts during initialization.",
     true);
  static SwitchOption interfaceInitCutsNo
    (interfaceInitCuts,
     "No",
     "Do not extract cuts during initialization.",
     false);


  interfaceScanCuts.rank(10.5);
  interfaceInitCuts.rank(10.6);

}

long MadGraphReader::numberOfEvents(string block) {
  long output(0);
  // Check for number of events in the file.
  string::size_type pos = block.find("##  Number of Events       :");
  if ( pos == string::npos )
    pos = block.find("#  Number of Events        :");
  if ( pos != string::npos ) {
    pos += 28;
    output = std::strtol(block.c_str() + pos, NULL, 0);
  }
  return output;
}
