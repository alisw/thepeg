#include "RandomHandler.h"

#include "ThePEG/Repository/UseRandom.h"


using namespace TheP8I;


	RandomHandler::RandomHandler(){
	}

	RandomHandler::~RandomHandler(){

	}
	
	RandomHandler::RandomHandler(bool walkerWeight) : _walkerWeight(walkerWeight) {

		// Recursion relations for addition of triplet and anti-triplet to a (p, q) multiplet
		addAntiTriplet.push_back(Plet(-1,0));
		addAntiTriplet.push_back(Plet(1,-1));
		addAntiTriplet.push_back(Plet(0,1));

		addTriplet.push_back(Plet(0,-1));
		addTriplet.push_back(Plet(-1,1));
		addTriplet.push_back(Plet(1,0));

	}

	double RandomHandler::KappaEnhancement(StringPipe& pipe){
	    // We loop over all pipes to calculate overlap
	    for(vector<StringPipe>::iterator sItr = _pipes.begin(); sItr!=_pipes.end(); ++sItr){
	    	// Only do something if we reach another pipe than this
	    	if(pipe._ymin == sItr->_ymin && pipe._ymax == sItr->_ymax && pipe._radius2 == sItr->_radius2){
	    		int m = floor(sItr->_internalOverlap.first + sItr->_externalOverlap.first + 0.5);
				int n = floor(sItr->_internalOverlap.second + sItr->_externalOverlap.second + 0.5);
			
				// This will become nan if one of the containers have a volume of zero.
				// In that case no overlap should be counted.	
				if(isnan(m+n)) return 1; 

				double p = 0;
				double q = 0;
				// These three cases are handled separately -- no overlap => no enhancement
				if((m==1 && n==0) || (m==0 && n==1) || (m==0 && n==0)){
					p = 1;
					q = 0;
				}
				else{
				// Do the full random walk
					// Start from this initial conf.
					Plet plet(0,0);
					int mstep = m;
					int nstep = n;
					// As long as there are still steps to be taken...
					while(mstep > 0 || nstep > 0){
							// Give weights to individual steps
							vector<double> we;
							we.clear();
							double prob = (double(plet.q()) + double(plet.p())) / (double(m) + double(n)); 
							if(_walkerWeight) we.push_back(1.0 - prob);
							else we.push_back(1.0);
							we.push_back(1.0);
							we.push_back(1.0);		
						// Take one at random if you can
						if(UseRandom::rnd() < 0.5 && mstep > 0){
							plet = plet.addRandomPlet(addTriplet,we,UseRandom::rnd());
							--mstep;
							}
						else if(nstep > 0){
							plet = plet.addRandomPlet(addAntiTriplet,we,UseRandom::rnd());
							--nstep;
							}
						}
				
					p = plet.p();
					q = plet.q();	
				}
				// TODO: Check why these configurations are sometimes returned
				if(isnan(p+q) || p + q == 0) return 1;

				// Should we throw the string?
				if(UseRandom::rnd() > (p + q)/(m + n)) return -999.0;
				//cout << p << " " << q << " " << m << " " << n << endl;
				double N = p + q;

				return (0.25*(N + 3 - p*q/N));
	    	}
	    }
	    cout << "Could not find pipe..." << endl;
	    AddPipe(pipe);
	    return KappaEnhancement(pipe);
	}

	void RandomHandler::SetEvent(vector<StringPipe> event){
	
		// When we set a new event, we will calculate all overlaps once and for all,
		// This will be stored as a std::vector of StringContainers _containers, which is a private
		// data member of the present class
		_pipes = event;
	    RecalculateOverlaps();

	}

	void RandomHandler::AddPipe(StringPipe& pipe){
		_pipes.push_back(pipe);
		RecalculateOverlaps();

	}

	void RandomHandler::RecalculateOverlaps(){
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

	bool RandomHandler::RemovePipe(StringPipe& pipe){
	    for(vector<StringPipe>::iterator sItr = _pipes.begin(); sItr!=_pipes.end(); ++sItr){
	    	if(pipe._ymin == sItr->_ymin && pipe._ymax == sItr->_ymax && pipe._radius2 == sItr->_radius2){
	    		_pipes.erase(sItr);
	    		RecalculateOverlaps();
	    		return true;
	    	}
		}
		return false;
	}

	void RandomHandler::clear(){
		_pipes.clear();
	}