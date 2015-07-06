// -*- C++ -*-
#ifndef THEP8I_StringPipe_H
#define THEP8I_StringPipe_H

#include "TheP8I/Config/Pythia8Interface.h"
#include "TheP8EventShapes.h"

namespace TheP8I {

using namespace ThePEG;

class StringPipe {
public:
    StringPipe();

   ~StringPipe();

	StringPipe(ColourSinglet* singlet, Length r0, TheP8EventShapes * es = NULL);

    pair<double,double> ExternalOverlap(StringPipe& other);

    Area GetVolume();

    Area OverlapArea(StringPipe& other);

    double IntevalOverlap(double min1, double min2, double max1, double max2);

    double OverlapY(StringPipe& other);

    double MaxpTRapidity();
    Energy MeanpT();

    Energy MaxpT();

    ColourSinglet * GetSingletPtr();

    pair<Length, Length> _originb;
    double _ymin, _ymax;
    Area _radius2;
    pair<double,double> _internalOverlap;
    pair<double,double> _externalOverlap;

    TheP8EventShapes * _es;
    ColourSinglet * theSinglet;

};

}
#endif