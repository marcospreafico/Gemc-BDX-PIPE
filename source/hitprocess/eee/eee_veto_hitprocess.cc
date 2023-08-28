// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "eee_veto_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


map<string, double> eee_veto_HitProcess::integrateDgt(MHit* aHit, int hitn) {
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int sector = identity[0].id;
	int veto_id = identity[1].id;
	int channel = identity[2].id;

	// Digitization Parameters

// initialize ADC and TDC

	int ADC1 = 0;
	int ADC2 = 0;
	int ADC3 = 0;
	int ADC4 = 0;
	int ADC5 = 0;
	int ADC6 = 0;
	int ADC7 = 0;
	int ADC8 = 0;
	int TDC1 = 4096;
	int TDC2 = 4096;
	int TDC3 = 4096;
	int TDC4 = 4096;
	int TDC5 = 4096;
	int TDC6 = 4096;
	int TDC7 = 4096;
	int TDC8 = 4096;


	///EEE digitization
	if (veto_id == 1000 || veto_id == 1001 || veto_id == 1002) {
		double sSizeX = aHit->GetDetector().dimensions[0];
		double sSizeY = aHit->GetDetector().dimensions[1];
		double sSizeZ = aHit->GetDetector().dimensions[2];

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4ThreeVector> Gpos = aHit->GetPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();

		unsigned int nsteps = Edep.size();
		double Etot = 0;

		for (unsigned int s = 0; s < nsteps; s++)
			Etot = Etot + Edep[s];

		double vX = 0.;
		double vXL = 0.;
		double vY = 0.;
		double vZ = 0.;
		double Thit = 0.;
		double sigmaX = 8.4; // 8.4 spread in mm USED 9.2mm as from JINST_044P
		double sigmaZ = 8.4; // 8.4 spread in mm  USED 15mm as from JINST_044P
		double sigmaT = 0.238; // 0.075 spread in ns 0.075 USED 238ps as from JINST_044P
		double deltaX = 0.;
		double deltaZ = 0.;
		double deltaT = 0.;
		double StripMult = -1.;
		double speed = 15.8; // 15.8 as in rec program different from 11.24 from NIM A593 (2008) 263
		double TimeLeft = 0, TimeRight = 0;

		if (Etot > 0 && times[0] > 0) {
			for (unsigned int s = 0; s < nsteps; s++)
				vX = vX + Gpos[s].x();
			for (unsigned int s = 0; s < nsteps; s++)
				vXL = vXL + Lpos[s].x();
			for (unsigned int s = 0; s < nsteps; s++)
				vY = vY + Gpos[s].y();
			for (unsigned int s = 0; s < nsteps; s++)
				vZ = vZ + Gpos[s].z();
			for (unsigned int s = 0; s < nsteps; s++) {
				Thit = Thit + times[s];
				//cout <<  "Thits: " <<times[s] <<    " z: " <<Lpos[s].z() << endl ;
			}
			// hit pos in cm with gaussian spread
			deltaX = G4RandGauss::shoot(0.,sigmaX);
			deltaZ = G4RandGauss::shoot(0.,sigmaZ);
			deltaT = G4RandGauss::shoot(0.,sigmaT);
			vXL = vXL / nsteps;
			double vXLsmeared = vXL + 2 * deltaX;
			if (sector == 1) { // hit in a strip

				StripMult = 1.;
				if (abs(deltaX) > (sSizeX + 7 + vXL)) StripMult = StripMult + 1.; // strip gap = 7mm
				if (abs(deltaX) > (sSizeX + 7 - vXL)) StripMult = StripMult + 1.; // strip gap = 7mm
				if (abs(deltaX) > (3 * sSizeX + 2 * 7 + vXL)) StripMult = StripMult + 1.; // strip gap = 7mm
				if (abs(deltaX) > (3 * sSizeX + 2 * 7 - vXL)) StripMult = StripMult + 1.; // strip gap = 7mm

			}

			if (sector == 2) { // hit in a gap

				StripMult = 0.;
				if (abs(deltaX) > (sSizeX + vXL)) StripMult = StripMult + 1.; //
				if (abs(deltaX) > (sSizeX - vXL)) StripMult = StripMult + 1.; //
				if (abs(deltaX) > (sSizeX + 25 + vXL)) StripMult = StripMult + 1.; // strip = 25mm
				if (abs(deltaX) > (sSizeX + 25 - vXL)) StripMult = StripMult + 1.; // strip = 25mm
				if (abs(deltaX) > (3 * sSizeX + 25 + vXL)) StripMult = StripMult + 1.; //strip = 25mm
				if (abs(deltaX) > (3 * sSizeX + 25 - vXL)) StripMult = StripMult + 1.; //strip = 25mm

			}
			// cout <<  " xL: " << vXL <<  " XLsm: " << vXLsmeared << " mult: " << StripMult << " check strip: "<< (sSizeX+7+vXL)<< " check gap: "<< (sSizeX+vXL)<< endl ;

			vX = (vX / (nsteps) + deltaX) / 10.; // in cm
			vY = (vY / (nsteps)) / 10.; // in cm
			vZ = (vZ / (nsteps) + deltaZ) / 10.; // in cm
			Thit = (Thit / (nsteps) + deltaT); // in ns
			TimeLeft = (sSizeZ / 10. - vZ) / speed;
			TimeRight = (sSizeZ / 10. + vZ) / speed;
		}

		// mimiking the trigger window

		//cout <<  " x: " << vX <<  " y: " << vY << " z: " << vZ <<  " T: " << Thit <<  " TL: " << TimeLeft <<  " TR: " << TimeRight <<  " Strip: " << channel <<  endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
		ADC1 = vX * 10000.;	        //in um (ADCX are INT!)
		ADC2 = vY * 10000.;	        //in um
		ADC3 = vZ * 10000; //in um
		ADC4 = StripMult;
		TDC1 = Thit * 1000; //in ps time from geant
		// output compatible with rec program
		int outmacro = 1; //1 to have output compatibele with macro 0 to have out compatible with reconstruction
		if (outmacro == 0) {
			ADC1 = StripMult * 2; // Nhits in the event
			ADC3 = channel; // Strip Left
			ADC4 = channel + 24; // Strip Right
			//int RefTime=89190+7836-761; // derived by tracks hitting the center of the middle chamber 3368
			int RefTime = 7836 - 761; // derived by tracks hitting the center of the middle chamber 3368
			TDC1 = Thit * 1000 - RefTime; //in ps time from geant
			TDC3 = (Thit + TimeLeft) * 1000 - RefTime; // hit time propagated to the L side
			TDC4 = (Thit + TimeRight) * 1000 - RefTime; // hit time propagated to the R side
			//cout <<  " x: " << vX <<  " y: " << vY << " z: " << vZ <<  " TDC: " << TDC1  <<  " Nhit: " << ADC1<<  " TL: " << TDC3 <<  " TR: " << TDC4 <<  " StripL: " << ADC3<<  " StripR: " << ADC4 <<  endl ;
		} //end out compatible with reconstruction
	}

	dgtz["hitn"] = hitn;
	dgtz["sector"] = sector;
	dgtz["veto"] = veto_id;
	dgtz["channel"] = channel;
	dgtz["adc1"] = ADC1;        // output in pe
	dgtz["adc2"] = ADC2;        //deposited energy in keV
	dgtz["adc3"] = ADC3;        // ignore
	dgtz["adc4"] = ADC4;        // ignore
	dgtz["adc5"] = ADC5;        // ignore
	dgtz["adc6"] = ADC6;        // ignore
	dgtz["adc7"] = ADC7;        // ignore
	dgtz["adc8"] = ADC8;        // ignore

	dgtz["tdc1"] = TDC1;        // output in ps
	dgtz["tdc2"] = TDC2;        // ignore
	dgtz["tdc3"] = TDC3;        // ignore
	dgtz["tdc4"] = TDC4;        // ignore
	dgtz["tdc5"] = TDC5;        // ignore
	dgtz["tdc6"] = TDC6;        // ignore
	dgtz["tdc7"] = TDC7;        // ignore
	dgtz["tdc8"] = TDC8;        // ignore

	return dgtz;
}

vector<identifier> eee_veto_HitProcess::processID(vector<identifier> id, G4Step *step, detector Detector) {
	id[id.size() - 1].id_sharing = 1;
	return id;
}



double eee_veto_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks) {
//Example of Birk attenuation law in organic scintillators.
//adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
//
// Taken from GEANT4 examples advanced/amsEcal and extended/electromagnetic/TestEm3
//
	double response = destep;
	if (birks * destep * stepl * charge != 0.) {
		response = destep / (1. + birks * destep / stepl);
	}
	return response;
}

double eee_veto_HitProcess::BirksAttenuation2(double destep, double stepl, int charge, double birks) {
//Extension of Birk attenuation law proposed by Chou
// see G.V. O'Rielly et al. Nucl. Instr and Meth A368(1996)745
//
//
	double C = 9.59 * 1E-4 * mm * mm / MeV / MeV;
	double response = destep;
	if (birks * destep * stepl * charge != 0.) {
		response = destep / (1. + birks * destep / stepl + C * pow(destep / stepl, 2.));
	}
	return response;
}

map<string, vector<int> > eee_veto_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map<string, vector<int> > MH;

	return MH;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> eee_veto_HitProcess::electronicNoise() {
	vector<MHit*> noiseHits;

// first, identify the cells that would have electronic noise
// then instantiate hit with energy E, time T, identifier IDF:
//
// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

// push to noiseHits collection:
// noiseHits.push_back(thisNoiseHit)

	return noiseHits;
}

// - charge: returns charge/time digitized information / step
map<int, vector<double> > eee_veto_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map<int, vector<double> > CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double eee_veto_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

