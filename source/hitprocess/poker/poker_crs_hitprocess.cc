// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "poker_crs_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> poker_crs_HitProcess::integrateDgt(MHit* aHit, int hitn) {
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int sector = identity[0].id;
	int xch = identity[1].id;
	int ych = identity[2].id;
        int zch = identity[3].id;

	
	// sensor parameters
    double integration_frac =1;
	double sensor_surface_crs = pow(0.6 * cm, 2);
 	double sensor_qe_crs = 0.22; // consider only 25um sipm
    double sensor_pixel_area = pow(0.0025*cm,2);
    double sensor_n_pixel = sensor_surface_crs/sensor_pixel_area;
 	double optical_coupling = 0.9;

 	// PbWO4 parameters
    double att_length_crs = 60000 * cm; // compatible with NO ATT Lenght as measured for cosmic muons
    double light_yield_crs = 310 * integration_frac * (1. / MeV);
	double length_crs; //(in mm)
	double sside_crs;
	double lside_crs;
    double redout_surface_crs;
	length_crs = 2 * aHit->GetDetector().dimensions[2];
	sside_crs = 2 * aHit->GetDetector().dimensions[0];
	lside_crs = 2 * aHit->GetDetector().dimensions[1];
	redout_surface_crs = sside_crs * lside_crs * mm * mm;

	double light_coll_crs = sensor_surface_crs / redout_surface_crs;
	if (light_coll_crs > 1) light_coll_crs = 1.;
	double etotL_crs = 0; //L= Large side redout

	double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();

	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> Dx = aHit->GetDx();


	vector<int> charge = aHit->GetCharges();
	vector<G4double> times = aHit->GetTime();
	//vector<string> theseMats = aHit->GetMaterials();

	unsigned int nsteps = Edep.size();
	double Etot_crs = 0;

	double peL_int_crs;
	double peL_crs = 0.;

	for (unsigned int s = 0; s < nsteps; s++) {
		Etot_crs = Etot_crs + Edep[s];
	}

	double Etot_B_crs = 0;
	double Etot_noB_crs = 0.;
	double peL_crs_poisson =0.;
	double peL_crs_sat =0;
	double peL_crs_poisson_sat =0;
	if (Etot_crs > 0) {
		for (unsigned int s = 0; s < nsteps; s++) {

		  //double dLeft_crs = length_crs / 2 + Lpos[s].z();            //Downstream (SIPM position )
			double Edep_B_crs = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
			Etot_B_crs = Etot_B_crs + Edep_B_crs;
			etotL_crs = etotL_crs + Edep_B_crs / 2 * exp(-length_crs / att_length_crs);

		}

		// input photon signal

		double n_photon = etotL_crs * light_yield_crs *light_coll_crs*optical_coupling;

		// number of pe
		peL_crs = n_photon* sensor_qe_crs;
        peL_crs_poisson = G4Poisson(peL_crs);

    	// number of pe - saturation effect;

    	peL_crs_sat = sensor_n_pixel*(1-exp(-n_photon*sensor_qe_crs/sensor_n_pixel));
    	peL_crs_poisson_sat = G4Poisson(peL_crs_sat);
	//

	}
	// closes (Etot > 0) loop
	if (verbosity > 4) {
	  cout << log_msg << " xch: " << xch << ", ych: " << ych <<", zch: "<<zch<<endl;
		cout << log_msg << " Etot=" << Etot_crs / MeV << endl;
	}

	dgtz["hitn"] = hitn;
	dgtz["sector"] = sector;
	dgtz["xch"] = xch;
	dgtz["ych"] = ych;
    dgtz["zch"] = zch;

	dgtz["Npe"] = peL_crs_poisson;	  // number of pe
	dgtz["Npe_sat"] = peL_crs_poisson_sat; // number of pe - saturation effect
	dgtz["Etot_B"] = Etot_B_crs;  // deposited energy with Birks
	dgtz["Etot_noB"] =  Etot_crs;



	return dgtz;
}

vector<identifier> poker_crs_HitProcess::processID(vector<identifier> id, G4Step *step, detector Detector) {
	id[id.size() - 1].id_sharing = 1;
	return id;
}

double poker_crs_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks) {
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

double poker_crs_HitProcess::BirksAttenuation2(double destep, double stepl, int charge, double birks) {
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

map<string, vector<int> > poker_crs_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map<string, vector<int> > MH;

	return MH;
}

//double crs_HitProcess::WaveForm(double npe, double time)
double* poker_crs_HitProcess::WaveForm(double npe, double* time) {
	double c = exp(-2.);
//    double Time;
	double t; // time in usec
	double WF;
	double y;
	double rr;
	int it;
	int Nch_digi = 800; //Number of cjannel for the digitizer
	static double WFsample[1000]; //Needs to be >  Nch_digi+size of the response to the single pe
	double smp_t = 4. / 1000.; // Assuming fADC sampling at 250 MHz 1sample every 4ns

	// double p[6] = {0.14,-3.5,2.5,-2.,0.5,-1.2};
	double p[6] = { 0., 0.680, 0.64, 3.34, 0.36, 0. }; // Babar CsI paprameters: fast component(in us), % fast, slow comp(in us), % slow
	// double p1[6] = {0.33,-0.04,3.45,-0.05,2.5,-0.045};

	double tau = 15.; // ampli response time constant (in ns)
	double t0 = 0.01; // t0 starting time (in ns)
	double area = (tau / c / 2.);
	double A = 1. / area; // amplitude at mnax (55.41 to have it normalized to integral=1, otherwise the max is at 1)
//    double threshold=10.*1./area/smp_t/1000.; //time threshold in pe - 1/55.41/smp_t*1000. is the funct max -

	double t_spread = 1. * 0.000; // pream time spread in us
	double A_spread = 1. * 0.05 * A; // pream amp spread (in fraction of 1pe amplitude = A)
	double func = 0.;
	// Building the waveform
	for (unsigned int s = 0; s < 1000; s++) {
		WFsample[s] = 0;
	}
	// Building the response to a single pe (preamps response)
	static double AmpWF[80];
	static int isFirst = 1;
	static double frac;// fraction of pe in Nch_digi

	if (isFirst) {
		for (unsigned int s = 0; s < 80; s++) {
			t = 1000. * s * smp_t;
			// parametrization of preamp out time is in ns (rise ~10ns decay~80ns) sampled in 160ns or 40 samples
			//func=1./411.5*((1-exp(p1[0]+p1[1]*t))*exp(p1[2]+p1[3]*t)+exp(p1[4]+p1[5]*t)));
			func = (t - t0) * (t - t0) * exp(-(t - t0) / tau) * A / (4 * tau * tau * c) * 0.5 * (abs(t - t0) / (t - t0) + 1);
			// spreading amplitude by apli noise
			AmpWF[s] = smp_t * 1000. * func;
		}
		frac = 1 - ((p[2] * exp(-smp_t * Nch_digi / p[1]) + p[4] * exp(-smp_t * Nch_digi / p[3])));// fraction of pe in Nch_digi
		isFirst = 0;
	}

	// fraction of pe in Nch_digi

	int Npe = frac * npe;

	for (unsigned int s = 1; s <= Npe; s++) {
		y = 1.;
		WF = 0.;
		while (y > WF) {
			rr = (rand() % 1000000 + 1) / 1000000.; // rnd number between 0-1
			t = Nch_digi * smp_t * rr; // extracting over 5000 samples range (5000x4ns=20us)
			//WF= 1./5.15*((1-exp(p[0]+p[1]*t))*exp(p[2]+p[3]*t)+exp(p[4]+p[5]*t));
			WF = (p[2] / p[1] * exp(-t / p[1]) + p[4] / p[3] * exp(-t / p[3])) / (p[2] / p[1] + p[4] / p[3]);
			rr = (rand() % 10000000 + 1) / 10000000.; // rnd number between 0-1
			y = rr;
			//  cout << "WF " << WF   << " rnd " << y  << endl;
		}
		// spreading time and amplitude of the ampli signal
		t = G4RandGauss::shoot(t, t_spread);
		if (t < 0.) t = 0.;
		it = t / smp_t;

		for (unsigned int s = 0; s < 80; s++) {
			t = 1000. * s * smp_t;
			func = AmpWF[s];
			func = G4RandGauss::shoot(func, A_spread);
			if ((s + it) < Nch_digi) WFsample[s + it] = WFsample[s + it] + func;
		}
	}

	// mimicking a CF discriminatorm at 1/3 of the max signal
	*time = 0.;
	double time_max = -100;
	int s = 0;
	int s_time_max = 0;
	while (time_max < WFsample[s]) {
		time_max = 1 / 2. * (WFsample[s + 1] + WFsample[s]);
		s_time_max = s;
		*time = 1000. * smp_t * s_time_max / 3.;
		s++;
	}
	// cout<<s_time_max<<"  "<< time_max<< "  "<<*time <<endl;

	/* // mimicking a FixedT discriminatorm
	 for(unsigned int s=0; s<1000; s++)
	 {
	 //cout << s  << " " <<  WFsample[s] << endl ;
	 //look for the max

	 if(WFsample[s]>threshold)
	 {*time=1000.*s*smp_t; //time in ns
	 break;
	 }
	 }
	 */
	return WFsample;

}

double* poker_crs_HitProcess::WaveFormPbwo(double npe, double* time_pbwo) {
	double c = exp(-2.);
	//    double Time;
	double t; // time in usec
	double WF;
	double y;
	double rr;
	int it;
	int Nch_digi = 800; //Number of channel for the digitizer
	static double WFsample[1000]; //Needs to be >  Nch_digi+size of the response to the single pe

	static int isFirst = 1;

	double smp_t = 4. / 1000.; // Assuming fADC sampling at 250 MHz 1sample every 4ns

	// double p[6] = {0.14,-3.5,2.5,-2.,0.5,-1.2};
	double p[6] = { 0., 0.00680, 0.64, 0.0334, 0.36, 0. }; // PbWO: fast component(in us), % fast, slow comp(in us), % slow
	// double p1[6] = {0.33,-0.04,3.45,-0.05,2.5,-0.045};

	double tau = 15.; // ampli response time constant (in ns)
	double t0 = 0.01; // t0 starting time (in ns)
	double area = (tau / c / 2.);
	double A = 1. / area; // amplitude at mnax (55.41 to have it normalized to integral=1, otherwise the max is at 1)
	//    double threshold=10.*1./area/smp_t/1000.; //time threshold in pe - 1/55.41/smp_t*1000. is the funct max -

	double t_spread = 1. * 0.000; // pream time spread in us
	double A_spread = 1. * 0.4 * A; // pream amp spread (in fraction of 1pe amplitude = A)
	double func = 0.;
	static double frac;	// fraction of pe in Nch_digi
	// Building the waveform
	for (unsigned int s = 0; s < 1000; s++) {
		WFsample[s] = 0;
	}
	// Building the response to a single pe (preamps response)
	static double AmpWF[80];
	if (isFirst) {
		for (unsigned int s = 0; s < 80; s++) {
			t = 1000. * s * smp_t;
			// parametrization of preamp out time is in ns (rise ~10ns decay~80ns) sampled in 160ns or 40 samples
			//func=1./411.5*((1-exp(p1[0]+p1[1]*t))*exp(p1[2]+p1[3]*t)+exp(p1[4]+p1[5]*t)));
			func = (t - t0) * (t - t0) * exp(-(t - t0) / tau) * A / (4 * tau * tau * c) * 0.5 * (abs(t - t0) / (t - t0) + 1);
			// spreading amplitude by apli noise
			AmpWF[s] = smp_t * 1000. * func;
		}
		frac = 1 - ((p[2] * exp(-smp_t * Nch_digi / p[1]) + p[4] * exp(-smp_t * Nch_digi / p[3])));	// fraction of pe in Nch_digi
		isFirst = 0;
	}

	//int mNpe = int(frac * npe);
	int mNpe = G4Poisson(frac * npe);
	for (unsigned int s = 1; s <= mNpe; s++) {
		y = 1.;
		WF = 0.;
		while (y > WF) {
			rr = (rand() % 1000000 + 1) / 1000000.; // rnd number between 0-1
			t = Nch_digi * smp_t * rr; // extracting over 5000 samples range (5000x4ns=20us)
			//WF= 1./5.15*((1-exp(p[0]+p[1]*t))*exp(p[2]+p[3]*t)+exp(p[4]+p[5]*t));
			WF = (p[2] / p[1] * exp(-t / p[1]) + p[4] / p[3] * exp(-t / p[3])) / (p[2] / p[1] + p[4] / p[3]);
			rr = (rand() % 10000000 + 1) / 10000000.; // rnd number between 0-1
			y = rr;
		}
		t = G4RandGauss::shoot(t, t_spread);
		if (t < 0.) t = 0.;
		it = t / smp_t;
		for (unsigned int s = 0; s < 80; s++) {
			t = 1000. * s * smp_t;
			func = AmpWF[s];
			func = G4RandGauss::shoot(func, A_spread);
			if ((s + it) < Nch_digi) WFsample[s + it] = WFsample[s + it] + func;
		}
	}

	// mimicking a CF discriminatorm at 1/3 of the max signal
	*time_pbwo = 0.;
	double time_max = -100;
	int s = 0;
	int s_time_max = 0;
	while (time_max < WFsample[s]) {
		time_max = 1 / 2. * (WFsample[s + 1] + WFsample[s]);
		s_time_max = s;
		*time_pbwo = 1000. * smp_t * s_time_max / 3.;
		s++;
	}
	// cout<<s_time_max<<"  "<< time_max<< "  "<<*time <<endl;

	/* // mimicking a FixedT discriminatorm
	 for(unsigned int s=0; s<1000; s++)
	 {
	 //cout << s  << " " <<  WFsample[s] << endl ;
	 //look for the max

	 if(WFsample[s]>threshold)
	 {*time=1000.*s*smp_t; //time in ns
	 break;
	 }
	 }
	 */
	return WFsample;

}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> poker_crs_HitProcess::electronicNoise() {
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
map<int, vector<double> > poker_crs_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map<int, vector<double> > CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double poker_crs_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

