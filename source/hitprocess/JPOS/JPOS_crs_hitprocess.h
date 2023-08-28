#ifndef JPOS_JPOS_crs_HitProcess_H
#define JPOS_JPOS_crs_HitProcess_H 1

// gemc headers
#include "HitProcess.h"

// Class definition
class JPOS_crs_HitProcess : public HitProcess
{
public:

	~JPOS_crs_HitProcess(){;}

	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);

	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);

	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);

	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);
	
	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);

	// creates the HitProcess
	static HitProcess *createHitClass() {return new JPOS_crs_HitProcess;}


	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();

};

#endif
