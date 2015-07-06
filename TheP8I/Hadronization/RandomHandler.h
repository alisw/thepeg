// -*- C++ -*-

#include <iostream>
#include <vector>
#include <math.h> 
#include "ThePEG/Repository/RandomGenerator.h"

using namespace std;
#ifndef THEP8I_Plet_H
#define THEP8I_Plet_H

class Plet{
public:
	Plet(int p, int q) : _p(p), _q(q) {
		// Construct SU(3) multiplet given quantum numbers p and q
		_disallowed.clear();
	}

	Plet add(int p, int q) {
		return(Plet(_p + p,_q + q));
	}

	Plet add(Plet t) {
		return add(t.p(),t.q());
	}

	double multiplicity(){
		// The multiplicity of an SU(3) multiplet
		if(_p<0 || _q<0) return 0;
		// Is it disallowed by hand?
		for(vector<Plet>::iterator it = _disallowed.begin(); it!=_disallowed.end(); ++it) 
			if(_p==it->p() && _q==it->q()) return 0;
		return 0.5*(_p + 1)*(_q + 1)*(_p + _q + 2);
	}

	Plet addRandomPlet(vector<Plet> plets, vector<double> weights, double ran){
	// Take a Plet at random from the plets vector, and add it to *this plet with corresponding weight
	// using multiplicities from above to calculate probability 	
	if(plets.size()!=weights.size()){
		cout << "Not equal number of weights and plets!" << endl;
		return *this;
	}

	vector<Plet> ret;
	for(vector<Plet>::iterator it = plets.begin(); it!=plets.end(); ++it) ret.push_back(add(*it)); 
	
	// Calculate normalization constant	
	double norm = 0;
	vector<double>::iterator it2 = weights.begin();
	for(vector<Plet>::iterator it = ret.begin(); it!=ret.end();++it, ++it2) norm+=it->multiplicity()*(*it2);

	// Sanity check	
	if(ret.size() == 0){
		cout << "Could not walk!" << endl;
		return *this;
	}

	// Calculate cumulated probability to walk to a given state, and check whether the random number is less than that.
	// If yes, return said state
	double cumprob = 0;
	it2 = weights.begin();
	for(vector<Plet>::iterator it = ret.begin(); it!=ret.end(); ++it, ++it2){
		cumprob += it->multiplicity()*(*it2)/norm;
			if(ran < cumprob) return *it;
	}

	cout << "We should never reach this point!" << endl;
	return ret.back();
}
	void disablePlet(Plet set){
		// Set disallowed plets, eg. singlets
		_disallowed.push_back(set);
	}

	int p(){
		return _p;
	}

	int q(){
		return _q;
	}

	int sum(){
		return _p+_q;
	}

private:
	int _p,_q;
	vector<Plet> _disallowed;
};

#endif


#ifndef THEP8I_RandomHandler_H
#define THEP8I_RandomHandler_H

#include "TheP8EventShapes.h"
#include "StringPipe.h"

namespace TheP8I {

using namespace ThePEG;

class RandomHandler {
public:
	RandomHandler();

	~RandomHandler();
	
	RandomHandler(bool walkerWeight);

	double KappaEnhancement(StringPipe& pipe);

	void SetEvent(vector<StringPipe> event);

	void AddPipe(StringPipe& pipe);

	void RecalculateOverlaps();

	bool RemovePipe(StringPipe& pipe);

	void clear();


private:
	vector<StringPipe> _pipes;
	vector<Plet> addAntiTriplet;
	vector<Plet> addTriplet;
	bool _walkerWeight;

};

}
#endif