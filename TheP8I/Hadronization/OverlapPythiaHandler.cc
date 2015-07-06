// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OverlapPythiaHandler class.
//

#include "OverlapPythiaHandler.h"

using namespace TheP8I;

OverlapPythiaHandler::OverlapPythiaHandler(){

}

OverlapPythiaHandler::~OverlapPythiaHandler(){
  for ( vector<PythiaPtr>::iterator itr = _overlapPythias.begin();
	itr != _overlapPythias.end(); ++itr )
    itr->nullify();
  _overlapPythias.clear();
}

OverlapPythiaHandler::OverlapPythiaHandler(HadronizationHandler * sf, vector<string>arguments) : _sf(sf){
  // Set the arguments
  for(vector<string>::iterator itr = arguments.begin(); itr!=arguments.end();++itr){
    if(itr->find("StringPT:sigma")!=string::npos){
      sigma = PythiaParameter(*itr);
    }

    else if(itr->find("StringZ:aLund")!=string::npos){
      a = PythiaParameter(*itr);
    }

    else if(itr->find("StringZ:bLund")!=string::npos){
      b = PythiaParameter(*itr);
    }

    else if(itr->find("StringFlav:probStoUD")!=string::npos){
      rho = PythiaParameter(*itr);
    }

    else if(itr->find("StringFlav:probSQtoQQ")!=string::npos){
      x = PythiaParameter(*itr);
    }
    
    else if(itr->find("StringFlav:probQQ1toQQ0")!=string::npos){
      y = PythiaParameter(*itr);
    }

    else if(itr->find("StringFlav:probQQtoQ")!=string::npos){
      xi = PythiaParameter(*itr);
    }

    else if(itr->find("OverlapStrings:fragMass")!=string::npos){
      double m = PythiaParameter(*itr);
      m2 = m*m;
    }

    else if(itr->find("OverlapStrings:baryonSuppression")!=string::npos){
      _bsparameter = PythiaParameter(*itr);
    }

    else{
      _pythiaSettings.push_back(*itr);
    }
  }
  // Sanity check of parameters
      if (rho < 0 || rho > 1 || x < 0 || x > 1 || y < 0 || y > 1 || xi < 0 ||
       xi > 1 || sigma < 0 || sigma > 1 || a < 0.0 || a > 2.0 || b < 0.2 || b > 2.0){
        cout << "Did you set up sensible initial Pythia 8 values? Remember:" << endl;
        cout << "0 < a < 2; 0.2 < b < 2; 0 < sigma < 1; 0 < xi < 1; 0 < y < 1; 0 < x < 1; 0 < rho < 1" << endl;
      }
}

void OverlapPythiaHandler::CreateSinglePythia(double h){

  if(!CalculateEffectiveParameters(h)) cout << "Unexpected error setting up parameters for overlap string mode!" << endl;
 /*  cout << "Settings and parameters for overlap string hadronization" << endl;
    cout << "bs = " << _bsparameter << " m_cut = " << sqrt(m2) << endl;
    cout << " h = " << h << " rho_eff = " << rho_eff << " xi_eff = " << xi_eff << " a_eff = " << a_eff << " b_eff = " << b_eff << " x_eff = " << x_eff << " y_eff = " << y_eff << endl;
  */
  vector<string> newPythiaSettings = _pythiaSettings;
  newPythiaSettings.push_back("StringPT:sigma = " + convert(sigma_eff));
  newPythiaSettings.push_back("StringZ:aLund = " + convert(a_eff));
  newPythiaSettings.push_back("StringZ:bLund = " + convert(b_eff));
  newPythiaSettings.push_back("StringFlav:probStoUD = " + convert(rho_eff));
  newPythiaSettings.push_back("StringFlav:probSQtoQQ = " + convert(x_eff));
  newPythiaSettings.push_back("StringFlav:probQQ1toQQ0 = " +  convert(y_eff));
  newPythiaSettings.push_back("StringFlav:probQQtoQ = " + convert(xi_eff));  
  
  PythiaPtr tmp;
  _overlapPythias.push_back(tmp);
  _overlapPythias.back().pPtr = new Pythia8Interface();
  _overlapPythias.back().pPtr->init(*_sf,newPythiaSettings);
  _overlapPythias.back().h = h;
  ComparePythias comparePythias;
  sort(_overlapPythias.begin(),_overlapPythias.end(),comparePythias);
}


Pythia8Interface * OverlapPythiaHandler::GetPythiaPtr(double h, bool recycle){
	// Broken calculations could give h = nan. This recursive function would thus 
	// loop infinitely, as nothing == nan. In that case, throw a Veto, and discard the event
	// after making sufficient noise.
	if(isnan(h)){
		cout << "OverlapPythiaHandler::GetPythiaPtr: h is nan. Fix preceeding calculation before running again." << endl;
	    throw Veto();
	}
	// Option to delete all unused Pythia ptrs every ef. 10000 events.
	// Usefull for Heavy Ion, which contains rare circumstances
  if(recycle){
      vector<PythiaPtr> newVec;
      for(vector<PythiaPtr>::iterator itr = _overlapPythias.begin(); itr!=_overlapPythias.end();++itr){
        if ( itr->h !=h ){
	  itr->nullify( );
        } 
        else{
          newVec.push_back(*itr);
        } 
      }
      _overlapPythias = newVec;
    }
   // Check if we already have  
  for(vector<PythiaPtr>::iterator itr = _overlapPythias.begin(); itr!=_overlapPythias.end();++itr){
    if(itr->h==h){
      itr->used();
      return itr->getpPtr();
    }
  }
  //cout << "Creating new Pythia with h = " << h << endl;
  CreateSinglePythia(h);
  return GetPythiaPtr(h);
}

 vector<double> OverlapPythiaHandler::GetPythiaParameters(double h){
   vector<double> ret;
   ret.clear();
   if(!CalculateEffectiveParameters(h)) cout << "Something went wrong calculating effective Pythia parameters!" << endl;
    ret.push_back(h); 
    ret.push_back(a_eff);   
    ret.push_back(b_eff);   
    ret.push_back(rho_eff);   
    ret.push_back(xi_eff);   
    ret.push_back(x_eff);   
    ret.push_back(y_eff);   
    return ret;
 }


bool OverlapPythiaHandler::CalculateEffectiveParameters(double h){
    if (h <= 0) return false;

    double hinv = 1.0 / h;

    rho_eff = pow(rho, hinv);
    x_eff = pow(x, hinv);
    y_eff = pow(3.0 * y, hinv) / 3.0;
    sigma_eff = sigma * sqrt(h);

    double C1 = _bsparameter * xi * (2 + rho) / (3 + 2*x*rho + 9*y + 6*x*rho*y + x*x*rho*rho +
          3*y*x*x*rho*rho);
    xi_eff = (pow(C1, hinv)/_bsparameter) * (3 + 2*x_eff*rho_eff + 9*y_eff +
            6*x_eff*rho_eff*y_eff + x_eff*x_eff*rho_eff*rho_eff +
            3*y_eff*x_eff*x_eff*rho_eff*rho_eff) / (2 + rho_eff);
    if (xi_eff > 1.0) xi_eff = 1.0;

    C1 = b / (2 + rho);
    b_eff = C1 * (2 + rho_eff);
    if (b_eff < 0.2) b_eff = 0.2;
    if (b_eff > 2.0) b_eff = 2.0;

    double N = IFragmentationF(a, b);
    double N_eff = IFragmentationF(a, b_eff);
    int s = sign(N - N_eff);
    double da = 0.1;
    a_eff = a - s * da;
    do {
      N_eff = IFragmentationF(a_eff, b_eff);
      if (sign(N - N_eff) != s) {
        s = sign(N - N_eff);
        da /= 10.0;
      }
      a_eff -= s * da;
      if (a_eff < 0.0) {a_eff = 0.0; break;}
      if (a_eff > 2.0) {a_eff = 2.0; break;}
    } while (da > 0.00001);

    return true;
}

double OverlapPythiaHandler::IFragmentationF(double a, double b){
    const int N = 100000;
    const double dz = 1.0 / N; //increment of z variable
    double Integral = 0.0;
    for (double z = dz; z < 1; z += dz) {
      //numerical integral calculation
      Integral += pow(1 - z, a) * exp(- b * m2 / z) / z;
    }
    return Integral * dz;
  }

double OverlapPythiaHandler::PythiaParameter(string par){
      size_t found = par.find("=");
      string number = par.substr(found+2,par.size());
      string::iterator end_pos = remove(par.begin(), par.end(), ' ');
       return atof(number.c_str());
}

int OverlapPythiaHandler::sign(double num){
    return (num < 0) ? -1 : 1;
  }
  
string OverlapPythiaHandler::convert(double d) {
  ostringstream os;
  os << d;
  return os.str();
}
