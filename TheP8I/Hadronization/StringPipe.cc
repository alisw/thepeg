#include "StringPipe.h"

using namespace TheP8I;

StringPipe::StringPipe() {

}

StringPipe::~StringPipe(){

}

StringPipe::StringPipe(ColourSinglet* singlet, Length r0, TheP8EventShapes * es) : _externalOverlap(make_pair(0,0)), _es(es), theSinglet(singlet) {
    		tcPPtr pmin = singlet->partons()[0];
    		tcPPtr pmax = singlet->partons()[0];
    		for (tcPVector::iterator pItr = singlet->partons().begin(); pItr!=singlet->partons().end();++pItr) {
    		  	// Find the lowest and the highest rapidities to span the full string between
    			double y = _es ? _es->yT((*pItr)->momentum()) : (*pItr)->rapidity();	

    			if(y<(_es ? _es->yT(pmin->momentum()) : pmin->rapidity())) pmin = *pItr;
    		  	else if(y>(_es ? _es->yT(pmax->momentum()) : pmax->rapidity())) pmax = *pItr;	
    		}
    		_originb.first = (pmin->vertex().x() +  pmax->vertex().x()) / 2;
    		_originb.second = (pmin->vertex().y() + pmax->vertex().y()) / 2;
    		_ymin = _es ? _es->yT(pmin->momentum()) : pmin->rapidity();
    		_ymax = _es ? _es->yT(pmax->momentum()) : pmax->rapidity();

    		pair<Area,Area> Vcs = make_pair(0*femtometer*femtometer,0*femtometer*femtometer);
    		tcPVector::iterator pminus = singlet->partons().begin();
            // Calculate pipe radius using distance from parton to mean line 
            // between partons with max and min probability
    		for (tcPVector::iterator pItr = singlet->partons().begin(); pItr!=singlet->partons().end();++pItr) {
    			// Grow the pipe radius
    			Length dx = (*pItr)->vertex().x() - _originb.first;
    			Length dy = (*pItr)->vertex().y() - _originb.second;
    			_radius2 += dx*dx + dy*dy; 

    			// Grow the volume of the colour singlet
    			if(pItr!=singlet->partons().begin()){
    				// Length dbx = (*pItr)->vertex().x() - (*pminus)->vertex().x();
    				// Length dby = (*pItr)->vertex().y() - (*pminus)->vertex().y();

    				double dy = _es ? _es->yT((*pItr)->momentum()) - _es->yT((*pminus)->momentum()) :
    						 (*pItr)->rapidity() - (*pminus)->rapidity();

    				// ORIGINAL BELOW -- THIS IS CHANGED PER 20.05.2014, pipe ends should now be parallel
                    // to impact parameter plane. Below they are allowed not to be.         
                       if(dy>0) Vcs.first += r0 * r0 * abs(dy);
                       if(dy<0) Vcs.second += r0 * r0 * abs(dy);     
                    // Area db2 = dbx * dbx + dby * dby;
    				// if(dy>0) Vcs.first += (r0 * r0 * dy * dy) / (sqrt(db2/(r0*r0) +dy*dy));
    				// if(dy<0) Vcs.second += (r0 * r0 * dy * dy) / (sqrt(db2/(r0*r0) +dy*dy));
    				++pminus;
    			}
 			}
            // Use radius squared instead of volume to save us from a factor of pi.
 			_radius2 /= (singlet->partons().end() - singlet->partons().begin());
 			_radius2  += r0*r0;

            // Overlap is given in both directions
 			_internalOverlap.first = Vcs.first / (_radius2 * (_ymax - _ymin));
 			_internalOverlap.second = Vcs.second / (_radius2 * (_ymax - _ymin));

    	}

pair<double,double> StringPipe::ExternalOverlap(StringPipe& other){
        if(this==&other || other.GetVolume() == 0*femtometer*femtometer){
            return make_pair(0,0);
       }
               return make_pair(other._internalOverlap.first * OverlapY(other) * OverlapArea(other) / other.GetVolume(),
                other._internalOverlap.second * OverlapY(other) * OverlapArea(other) / other.GetVolume());
    }

Area StringPipe::GetVolume(){
         return _radius2*M_PI*(_ymax - _ymin);
    }

Area StringPipe::OverlapArea(StringPipe& other){
        // This is the overlap between two cylinders. Purely geometrical. 
        Area d2 = (_originb.first -  other._originb.first) * (_originb.first -  other._originb.first) +
                 (_originb.second - other._originb.second) * (_originb.second - other._originb.second);
        Length d = sqrt(d2);
        Length r = sqrt(_radius2);
        Length R = sqrt(other._radius2);
        Qty<4,0,0,1,1,1> argu = (-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R);
        // Save us from dividing by  zero
        if(argu <= 0*femtometer*femtometer*femtometer*femtometer) return 0*femtometer*femtometer;
        return _radius2 * acos((d2 + _radius2 - other._radius2)/(2*d*r)) + other._radius2 * acos((d2 + other._radius2 - _radius2)/(2*d*R)) - 
                    0.5 * sqrt(argu);

    }

double StringPipe::IntevalOverlap(double min1, double min2, double max1, double max2){
            // Helper function to calculate overlap between two 1D intervals.
            // This could very well exist somewhere in ThePEG already.
             if(max1<min2||max2<min1||(min1==min2&&max1==max2))
                 return 0.;
            double ret = 0;
            if(min1<min2){
                 if(max1<max2) ret = abs(max1 - min2);
                 else {
                     ret = abs(max2 - min2);
                 }
             }
            else {
                if(max2<max1) ret = abs(max2 - min1);
                else {
                    ret = abs(max1 - min1);
                   }
            }
            return ret;
    }

double StringPipe::MaxpTRapidity(){
        // Gives the rapidity of the parton with the maximal pT
        Energy pT = 0*GeV;
        double y = 0;        
        for (tcPVector::iterator pItr = theSinglet->partons().begin(); pItr!=theSinglet->partons().end();++pItr) {
            if(pT < (*pItr)->momentum().perp()){
                pT = (*pItr)->momentum().perp();
                y = (*pItr)->rapidity();
            }
        }
        return y;

}    


Energy StringPipe::MeanpT(){
        Energy pT = 0*GeV;
        
        for (tcPVector::iterator pItr = theSinglet->partons().begin(); pItr!=theSinglet->partons().end();++pItr) {
            pT += (*pItr)->momentum().perp();
        }

        return pT/(distance(theSinglet->partons().begin(),theSinglet->partons().end()));
}

Energy StringPipe::MaxpT(){
        Energy pT = 0*GeV;
        
        for (tcPVector::iterator pItr = theSinglet->partons().begin(); pItr!=theSinglet->partons().end();++pItr) {
            if(pT < (*pItr)->momentum().perp()) pT = (*pItr)->momentum().perp();
        }
        return pT;
}

double StringPipe::OverlapY(StringPipe& other){
            return IntevalOverlap(_ymin, other._ymin, _ymax, other._ymax);
    }

ColourSinglet * StringPipe::GetSingletPtr(){
    return theSinglet;
}