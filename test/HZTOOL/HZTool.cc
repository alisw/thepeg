// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HZTool class.
//

#include "HZTool.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "HepMC/GenEvent.h"
#include "ThePEG/Vectors/HepMCConverter.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "HepMC/IO_HEPEVT.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HZTool.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace ThePEG {

template<>
struct HepMCTraits<HepMC::GenEvent>:
    public HepMCTraitsBase<HepMC::GenEvent,HepMC::GenParticle,
                           HepMC::GenVertex,HepMC::Polarization,
			   HepMC::PdfInfo> {};

}

extern "C" {

  void hztoolinit_(const char *, int);
  void hztoolsetxsec_(const double &, const double &);
  void hztoolfinish_();
  void hztoolredirect_(const char *, int);
  void hztoolfixdaughters_();
  void hztoolfixbeams_();
  void hztoolfixboson_(const double &, const double &, const double &,
		       const double &, const double &);
  void hztoolfixmirdir_();
  void hzfilhep_();
  void hz95108_(const int &);
  void hz96215_(const int &);
  void hz98076_(const int &);
  void hz98050_(const int &);
  void hz98051_(const int &);
  void hz98143_(const int &);
  void hzf01211e_(const int &);
  void hzf89201e_(const int &);
  void hzh9807014_(const int &);
  void hzh9807018_(const int &);
  void hzh9905024_(const int &);
  void hzh9907009_(const int &);
  void hzh9912022_(const int &);
  void hzh0001021_(const int &);
  void hzh0010026_(const int &);
  //  void hz0307080_(const int &);
  void hzh0412071_(const int &);

  void printhepevtp_();

}

using namespace ThePEG;

HZTool::HZAnaFn HZTool::getFunctionPointer(string name) {
  if ( name == "" ) return 0;
  else if ( name == "hz95108" ) return hz95108_;
  else if ( name == "hz96215" ) return hz96215_;
  else if ( name == "hz98076" ) return hz98076_;
  else if ( name == "hz98050" ) return hz98050_;
  else if ( name == "hz98051" ) return hz98051_;
  else if ( name == "hz98143" ) return hz98143_;
  else if ( name == "hzf01211e" ) return hzf01211e_;
  else if ( name == "hzf89201e" ) return hzf89201e_;
  else if ( name == "hzh9807014" ) return hzh9807014_;
  else if ( name == "hzh9807018" ) return hzh9807018_;
  else if ( name == "hzh9905024" ) return hzh9905024_;
  else if ( name == "hzh9907009" ) return hzh9907009_;
  else if ( name == "hzh9912022" ) return hzh9912022_;
  else if ( name == "hzh0001021" ) return hzh0001021_;
  else if ( name == "hzh0010026" ) return hzh0010026_;
  //  else if ( name == "hz0307080" ) return hz0307080_;
  else if ( name == "hzh0412071" ) return hzh0412071_;
  else return 0;
}


HZTool::~HZTool() {}

void HZTool::analyze(tEventPtr event, long ieve, int loop, int state) {
  sumweight += event->weight();
  HepMC::GenEvent * geneve =
    HepMCConverter<HepMC::GenEvent>::convert(*event, true);
  static HepMC::IO_HEPEVT converter;
  converter.set_trust_both_mothers_and_daughters(true);
//   static int count = 0;
//   if ( count++ < 10 ) {
//     cerr << *event;
//     geneve->print(cerr);
//     printhepevtp_();
//   }
  converter.write_event(geneve);
  addDISLepton(*event->primarySubProcess());
  hztoolfixbeams_();
  hzfilhep_();
  hztoolfixdaughters_();
  for ( int i = 0, N = functions.size(); i < N; ++i )
    if ( functions[i] ) (*functions[i])(2);
  delete geneve;
}

LorentzRotation HZTool::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void HZTool::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void HZTool::analyze(tPPtr) {}

void HZTool::addDISLepton(const SubProcess & sub) {
  tcPPtr in;
  tcPPtr out;
  if ( LeptonMatcher::Check(sub.incoming().first->data()) )
    in = sub.incoming().first;
  if ( LeptonMatcher::Check(sub.incoming().second->data()) ) {
    if ( in ) return;
    in = sub.incoming().second;
  }
  if ( !in ) return;
  long idin = in->id();
  for ( int i = 0, N = sub.outgoing().size(); i < N; ++i ) {
    long idout = sub.outgoing()[i]->id();
    if ( !LeptonMatcher::Check(idout) ) continue;
    if ( idin == idout ) out = sub.outgoing()[i];
  }
  if ( !out ) return;

  Lorentz5Momentum pgam = in->momentum() - out->momentum();
  hztoolfixboson_(pgam.x()/GeV, pgam.y()/GeV, pgam.z()/GeV,
		  pgam.t()/GeV, pgam.m()/GeV);

  if ( in->momentum().z() < 0*GeV )
    hztoolfixmirdir_();

}

void HZTool::dofinish() {
  AnalysisHandler::dofinish();
  hztoolsetxsec_(generator()->integratedXSec()/nanobarn, sumweight);
  for ( int i = 0, N = functions.size(); i < N; ++i )
    if ( functions[i] ) (*functions[i])(3);
  hztoolfinish_();
}

void HZTool::doinitrun() {
  sumweight = 0.0;
  AnalysisHandler::doinitrun();
  string file = ( filename().empty()? generator()->filename(): filename() )
    + ".hzlog";
  hztoolredirect_(file.c_str(), file.length());
  file = ( filename().empty()? generator()->filename(): filename() )
    + ".rz";
  hztoolinit_(file.c_str(), file.length());
  for ( int i = 0, N = functions.size(); i < N; ++i )
    if ( functions[i] ) (*functions[i])(1);
}

void HZTool::persistentOutput(PersistentOStream & os) const {
  os << theFilename << functionNames << sumweight;
}

void HZTool::persistentInput(PersistentIStream & is, int) {
  is >> theFilename >> functionNames >> sumweight;
  functions.resize(functionNames.size());
  for ( int i = 0, N = functionNames.size(); i < N; ++i ) {
    functions[i] = getFunctionPointer(functionNames[i]);
    if ( !functions[i] ) Throw<MissingFunction>()
      << "While reading '" << fullName() << "': The function '"
      << functionNames[i]
      << "' was not present in the installed HZTool version."
      << Exception::runerror;
  }
}

void HZTool::insertFunction(string name, int pos) {
  HZAnaFn fp = getFunctionPointer(name);
  if ( !fp ) Throw<MissingFunction>()
    << "The function '" << name
    << "' was not present in the installed HZTool version."
    << Exception::eventerror;
  pos = max(0, min(int(functionNames.size()), pos));
  functionNames.insert(functionNames.begin() + pos, name);
  functions.insert(functions.begin() + pos, fp);
}

void HZTool::setFunction(string name, int pos) {
  if ( pos < 0 || pos >= int(functionNames.size()) ) return;
  HZAnaFn fp = getFunctionPointer(name);
  if ( !fp ) Throw<MissingFunction>()
    << "The function '" << name
    << "' was not present in the installed HZTool version."
    << Exception::eventerror;
  functionNames[pos] = name;
  functions[pos] = fp;
}

void HZTool::delFunction(int pos) {
  if ( pos < 0 || pos >= int(functionNames.size()) ) return;
  functionNames.erase(functionNames.begin() + pos);
  functions.erase(functions.begin() + pos);
}

ClassDescription<HZTool> HZTool::initHZTool;
// Definition of the static class description member.

void HZTool::Init() {

  static ClassDocumentation<HZTool> documentation
    ("This class wraps the HZTool fortran library. The specified "
     "<interface>Functions</interface> of in the HZTool library will be "
     "called for each event (after the ThePEG::Event is first converted "
     "to a HepMC::GenEvent and then translated to the HEPEVT fortran common "
     "block). Note that only one HZTool AnalysisHandler can be used for "
     "a given EventGenerator, otherwise the result is undefined.");

  static Parameter<HZTool,string> interfaceFilename
    ("Filename",
     "The filename to which the HZTool histograms are written. If empty the "
     "name of the controlling EventGenerator is used instead. The standard "
     "'.rz' suffix is added to the name. Note that since the file is created "
     "within the cernlib hbook routines, the actual filename will be in lower "
     "case irrespectively of what is specified here. A file with the same "
     "name but with the suffix '.hzlog' will also be created containing all "
     "standard output from the HZTool library.",
     &HZTool::theFilename, "",
     true, false);

  static ParVector<HZTool,string> interfaceFunctions
    ("Functions",
     "The names of the HZTool analysis routines to be called in this analysis "
     "handler. Note that the routine names must be given in lower case "
     "letters.",
     &HZTool::functionNames, -1, "", "", "",
     true, false, Interface::nolimits,
     &HZTool::setFunction, &HZTool::insertFunction,
     &HZTool::delFunction, (vector<string>(HZTool::*)()const)(0),
     (string(HZTool::*)(int)const)(0), (string(HZTool::*)(int)const)(0),
     (string(HZTool::*)(int)const)(0), (vector<string>(HZTool::*)()const)(0));

}

