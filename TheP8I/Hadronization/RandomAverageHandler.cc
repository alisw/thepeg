#include "RandomAverageHandler.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace TheP8I;


	RandomAverageHandler::RandomAverageHandler(bool throwaway){

	}

	RandomAverageHandler::~RandomAverageHandler(){

	}


 	double RandomAverageHandler::KappaEnhancement(StringPipe& pipe){
	    for(vector<StringPipe>::iterator sItr = _pipes.begin(); sItr!=_pipes.end(); ++sItr){
	    	if(pipe._ymin == sItr->_ymin && pipe._ymax == sItr->_ymax && pipe._radius2 == sItr->_radius2){
	    		double m = floor(sItr->_internalOverlap.first + sItr->_externalOverlap.first + 0.5);
				double n = floor(sItr->_internalOverlap.second + sItr->_externalOverlap.second + 0.5);
	
				// Find out if the +1 should be left or right.
				// if(sItr->_internalOverlap.first > sItr->_internalOverlap.second) m+=1;
				// else n+=1;
			
				// This will become nan if one of the containers have a volume of zero.
				// In that case no overlap should be counted.	
				if(isnan(m+n)) return 1; 
				
				double p = 0;
				double q = 0;
				// These two^M^M^M three cases are handled separately
				if((m==1 && n==0) || (m==0 && n==1) || (m==0 && n==0)){
					p = 1;
					q = 0;
				}
				else{
					p = pow(m+n,17./24.)/2;
					q = p;
				}

				
				// Should we throw the string?
				if(UseRandom::rnd() > (p + q)/(m + n)) return -999.0;
				double N = p + q;

				return (0.25*(N + 3 - p*q/N));
	    	}
	    }
	    cout << "Could not find pipe..." << endl;
	    AddPipe(pipe);
	    return KappaEnhancement(pipe);
	}

	void RandomAverageHandler::SetEvent(vector<StringPipe> event){
	
	// When we set a new event, we will calculate all overlaps once and for all,
	// This will be stored as a std::vector of StringContainers _containers, which is a private
	// data member of the present class
	_pipes = event;
    RecalculateOverlaps();

	}

	void RandomAverageHandler::AddPipe(StringPipe& pipe){
		_pipes.push_back(pipe);
		RecalculateOverlaps();

	}

	void RandomAverageHandler::RecalculateOverlaps(){
	// Calculate all external overlaps to get the full number of overlapping strings per container
    for(vector<StringPipe>::iterator sItr = _pipes.begin(); sItr!=_pipes.end(); ++sItr){
    	sItr->_externalOverlap.first = 0;
    	sItr->_externalOverlap.second = 0;
		for(vector<StringPipe>::iterator sItr2 = _pipes.begin(); sItr2!=_pipes.end(); ++sItr2){
			sItr->_externalOverlap.first += (sItr->ExternalOverlap(*sItr2)).first;
			sItr->_externalOverlap.second += (sItr->ExternalOverlap(*sItr2)).second;
			}
		}

	}

	bool RandomAverageHandler::RemovePipe(StringPipe& pipe){
	    for(vector<StringPipe>::iterator sItr = _pipes.begin(); sItr!=_pipes.end(); ++sItr){
	    	if(pipe._ymin == sItr->_ymin && pipe._ymax == sItr->_ymax && pipe._radius2 == sItr->_radius2){
	    		_pipes.erase(sItr);
	    		RecalculateOverlaps();
	    		return true;
	    	}
		}
		return false;
	}

	void RandomAverageHandler::clear(){
		_pipes.clear();
	}