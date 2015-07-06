// -*- C++ -*-

#ifndef THEP8I_RandomAverageHandler_H
#define THEP8I_RandomAverageHandler_H

#include "TheP8EventShapes.h"
#include "StringPipe.h"

namespace TheP8I {

using namespace ThePEG;

class RandomAverageHandler {
public:
	RandomAverageHandler(bool throwaway);

	~RandomAverageHandler();
	
	double KappaEnhancement(StringPipe& pipe);

	void SetEvent(vector<StringPipe> event);

	void AddPipe(StringPipe& pipe);

	void RecalculateOverlaps();

	bool RemovePipe(StringPipe& pipe);

	void clear();

private:
	vector<StringPipe> _pipes;
};

}
#endif