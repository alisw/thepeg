// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ParameterHandler class.
//

#include "ParameterHandler.h"

using namespace TheP8I;

typedef map<string,double> PytPars;

ParameterHandler::ParameterHandler(){

}

ParameterHandler::~ParameterHandler(){

}

void ParameterHandler::init(double m2, double bsparameter, PytPars pars){
  _m2 = m2;
  _bsparameter = bsparameter;
  _parameters.clear();
  _parameters.insert(make_pair<double,PytPars>(1.0,pars));

  for (PytPars::iterator itr = pars.begin(); itr!=pars.end(); ++itr){

    if(itr->first.find("StringPT:sigma")!=string::npos) sigma = itr->second;

    else if(itr->first.find("StringZ:aLund")!=string::npos) a = itr->second;

    else if(itr->first.find("StringZ:bLund")!=string::npos) b = itr->second;

    else if(itr->first.find("StringFlav:probStoUD")!=string::npos) rho = itr->second;

    else if(itr->first.find("StringFlav:probSQtoQQ")!=string::npos) x = itr->second;
    
    else if(itr->first.find("StringFlav:probQQ1toQQ0")!=string::npos) y = itr->second;

    else if(itr->first.find("StringFlav:probQQtoQ")!=string::npos) xi = itr->second;

    else{
      cout << "Broken arrow!" << endl;
      }
    }
  // Sanity check of parameters
    if (rho < 0 || rho > 1 || x < 0 || x > 1 || y < 0 || y > 1 || xi < 0 ||
       xi > 1 || sigma < 0 || sigma > 1 || a < 0.0 || a > 2.0 || b < 0.2 || b > 2.0){
        cout << "Did you set up sensible initial Pythia 8 values? Remember:" << endl;
        cout << "0 < a < 2; 0.2 < b < 2; 0 < sigma < 1; 0 < xi < 1; 0 < y < 1; 0 < x < 1; 0 < rho < 1" << endl;
      }
}

PytPars ParameterHandler::GetEffectiveParameters(double h){  
   map<double,PytPars>::iterator it = _parameters.find(h);
   if(it!=_parameters.end()) return it->second;
   
   if(!CalculateEffectiveParameters(h)){
      cout << "Something went wrong calculating effective Pythia parameters!" << endl;
      return _parameters.find(1.0)->second;
    }
   if(!InsertEffectiveParameters(h)){
      cout << "Something went wrong inserting effective Pythia parameters!" << endl;
      return _parameters.find(1.0)->second;
    }
   return GetEffectiveParameters(h);
 }

bool ParameterHandler::InsertEffectiveParameters(double h){
    PytPars p;
    p.insert(make_pair<const string,double>("StringPT:sigma",sigma_eff));
    p.insert(make_pair<const string,double>("StringZ:aLund",a_eff));
    p.insert(make_pair<const string,double>("StringZ:bLund",b_eff));
    p.insert(make_pair<const string,double>("StringFlav:probStoUD",rho_eff));
    p.insert(make_pair<const string,double>("StringFlav:probSQtoQQ",x_eff));
    p.insert(make_pair<const string,double>("StringFlav:probQQ1toQQ0",y_eff));
    p.insert(make_pair<const string,double>("StringFlav:probQQtoQ",xi_eff));
    return (_parameters.insert(make_pair<double,PytPars>(h,p)).second);
}

bool ParameterHandler::CalculateEffectiveParameters(double h){
    if (h <= 0) return false;

    double hinv = 1.0 / h;

    rho_eff = pow(rho, hinv);
    x_eff = pow(x, hinv);
    y_eff = pow(y, hinv);
    sigma_eff = sigma * sqrt(h);

    double alpha = (1 + 2*x*rho + 9*y + 6*x*rho*y + 3*y*x*x*rho*rho)/(2 + rho);
    double alpha_eff = (1 + 2*x_eff*rho_eff + 9*y_eff + 6*x_eff*rho_eff*y_eff + 3*y_eff*x_eff*x_eff*rho_eff*rho_eff)/(2 + rho_eff);

    xi_eff = alpha_eff*_bsparameter*pow(xi/alpha/_bsparameter,hinv);

    if (xi_eff > 1.0) xi_eff = 1.0;

    b_eff = (2 + rho_eff)/(2 + rho) * b;
    if (b_eff < 0.2) b_eff = 0.2;
    if (b_eff > 2.0) b_eff = 2.0;

    double N = IFragmentationF(a, b);
    double N_eff = IFragmentationF(a, b_eff);
    int s = (N - N_eff) < 0 ? -1 : 1;
    double da = 0.1;
    a_eff = a - s * da;
    do {
      N_eff = IFragmentationF(a_eff, b_eff);
      if (((N - N_eff) < 0 ? -1 : 1) != s) {
        s = (N - N_eff) < 0 ? -1 : 1;
        da /= 10.0;
      }
      a_eff -= s * da;
      if (a_eff < 0.0) {a_eff = 0.0; break;}
      if (a_eff > 2.0) {a_eff = 2.0; break;}
    } while (da > 0.00001);

    return true;
}

double ParameterHandler::IFragmentationF(double a, double b){
    const int N = 100000;
    const double dz = 1.0 / N; //increment of z variable
    double Integral = 0.0;
    for (double z = dz; z < 1; z += dz) {
      //numerical integral calculation
      Integral += pow(1 - z, a) * exp(- b * _m2 / z) / z;
    }
    return Integral * dz;
  }

