// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "crs_hitprocess.h"

// Root utils
#include <TF1.h>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> crs_HitProcess::integrateDgt(MHit* aHit, int hitn) {
  
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
    
    int sector = identity[0].id;    // detector module
	int xch = identity[1].id;       // position of element in the module
	int ych = identity[2].id;
    int zch = identity[3].id;

    // Parameter for plastic
    double MeV2Pe = 60.; // taken from BDX-MINI data
    
    double time = 0.;
    double Etot = 0., Etot_B = 0.;
    
    double is_hit = 0.;  // This is the variable that tells us if there was a hit in the channel
    
    // TO BE IMPLEMENTED
    /*
     Ideally, we want to keep the electronics very simple to accomodate for the very large number of channels
     For this reason we can not digitalize fully the hit, but rather use some simpler evaluation of the charge
     The idea that Marco B. proposed is to use the same readout as for EPIC RICH, that has just a Time Over Threshold
     evaluation of the total deposited energy. For this reason the total energy has to be replaced with something similar
     */
    double TOT = 0.;
    
	
	// Get info about detector material to eveluate Birks effect
	double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();

	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> Dx = aHit->GetDx();

	// Charge for each step
	vector<int> charge = aHit->GetCharges();
	vector<G4double> times = aHit->GetTime();

	unsigned int nsteps = Edep.size();
        
	for (unsigned int s = 0; s < nsteps; s++) {
		Etot = Etot + Edep[s];
	}
    
	if (Etot > 0) {
		for (unsigned int s = 0; s < nsteps; s++) {
            
            double Edep_B = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
            
            Etot_B += Edep_B;
            
            double dLeft_crs = length_crs / 2 + Lpos[s].z();            //Downstream (SIPM position )
            double dRight_crs = length_crs / 2 - Lpos[s].z();            //Upstream
            
            // evaluate total energy and time for left readout
            etotL_crs += Edep_B / 2 * exp(-dLeft_crs / att_length_crs);
            timeL_crs += (times[s] + dLeft_crs / veff_crs) / nsteps;
            if (etotL_crs > 0.) {
                if (s == 0 || (time_min_crs[0] > (times[s] + dLeft_crs / veff_crs))) time_min_crs[0] = times[s] + dLeft_crs / veff_crs;
            }
            
            
            // evaluate total energy and time for right readout
            etotR_crs += Edep_B / 2 * exp(-dRight_crs / att_length_crs);
            timeR_crs += (times[s] + dRight_crs / veff_crs) / nsteps;
            if (etotR_crs > 0.) {
                if (s == 0 || (time_min_crs[1] > (times[s] + dRight_crs / veff_crs))) time_min_crs[1] = times[s] + dRight_crs / veff_crs;
            }
            
        }
        
		//      Right readout
		peR_crs = int(etotR_crs * light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs);
		peR_crs = G4Poisson(peR_crs);
        
        if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "CsI_Tl"){
            test = WaveForm(peR_crs, &tim);
        }else if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "G4_PbWO4"){
            test = WaveFormPbwo(peR_crs, &tim);
        }
		        
		for (unsigned int s = 0; s < Nsamp_int; s++) {
			peR_int_crs += test[s];
		}

        //      Left readout
        peL_crs = int(etotL_crs * light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs);
        peL_crs = G4Poisson(peL_crs);
        if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "G4_PbWO4"){
            test = WaveFormPbwo(peL_crs, &tim);
        } else if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "CsI_Tl"){
            test = WaveForm(peL_crs, &tim);
        }
        
        for (unsigned int s = 0; s < Nsamp_int; s++) {
            peL_int_crs = peL_int_crs + test[s];
        }
        
        // Save variables
        
        // Correct readout to be in energydouble sigfrac = 0;
        double  ts =  0.680, fs =  0.64, tl =  3.34, fl = 0.36; // fraction of long /short time; value of long/short time
        if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "G4_PbWO4"){
            ts = 0.00680; fs= 0.64; tl = 0.0334; fl =  0.36;
        } else if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "CsI_Tl"){
            ts =  0.680; fs =  0.64; tl =  3.34; fl = 0.36;
        }
        double digiframe =  Nsamp_int * 4. / 1000.;
        double sigfrac = 1 - (fs* exp(-digiframe / ts) + fl * exp(-digiframe / tl)); // fraction of signal contained in a digiframe digitalization window
        
        // energy measured = numbrer of phe / (light yield * attenuation * light was 2x before splitting * fraction measured in a fixed time window
        ADCR_crs = (peR_int_crs)/(light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs * 0.5 * sigfrac); // in MeV
        ADCL_crs = (peL_int_crs)/(light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs * 0.5 * sigfrac);
        
        //cout << peL_crs << " " << peL_int_crs/sigfrac << " ( was " << peL_int_crs << ")" << endl;
        
		TDCR_crs = int(tim) + ((time_min_crs[1] + T_offset_crs + G4RandGauss::shoot(0., sigmaTR_crs)) * tdc_conv_crs);
        TDCL_crs = int(tim) + ((time_min_crs[0] + T_offset_crs + G4RandGauss::shoot(0., sigmaTR_crs)) * tdc_conv_crs); // assigning to L the sipm2
        
		//Assigning to TDCB the usual timing seen by sipm1 (TDCR)
		TDCB = ((time_min_crs[1] + T_offset_crs + G4RandGauss::shoot(0., sigmaTR_crs)) * tdc_conv_crs);
        
        // overwrite right readout
        TDCR_crs = TDCL_crs;
        ADCR_crs = ADCL_crs;
    }
	// closes (Etot > 0) loop
    
	if (verbosity > 4) {
		cout << log_msg << " xch: " << xch << ", ych: " << ych;
		cout << log_msg << " Etot=" << Etot / MeV << endl;
		cout << log_msg << " TDCL=" << TDCL_crs << " TDCR=" << TDCR_crs << " ADCL=" << ADCL_crs << " ADCR=" << ADCR_crs << endl;
		//cout <<  log_msg << " TDCB=" << TDCB     << " TDCF=" << TDCF    << " ADCB=" << ADCB << " ADCF=" << ADCF << endl;
	}
    
	dgtz["hitn"] = hitn;
	dgtz["sector"] = sector;
	dgtz["xch"] = xch;
	dgtz["ych"] = ych;
    dgtz["zch"] = zch;
	dgtz["adcl"] = ADCL_crs;	  //
	dgtz["adcr"] = ADCR_crs;	  //SIPM 25um -> large size for matrix, small size for single
	dgtz["tdcl"] = TDCL_crs;	  //
	dgtz["tdcr"] = TDCR_crs;	  // as per ADCR_crs
	dgtz["adcb"] = Etot_B * 1000;  // deposited energy with Birks
	dgtz["adcf"] = Etot * 1000;
	dgtz["tdcb"] = TDCB * 1000.;	  //original time in ps
	dgtz["tdcf"] = 0;
	return dgtz;
}

vector<identifier> crs_HitProcess::processID(vector<identifier> id, G4Step *step, detector Detector) {
	id[id.size() - 1].id_sharing = 1;
	return id;
}

double crs_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks) {
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

double crs_HitProcess::BirksAttenuation2(double destep, double stepl, int charge, double birks) {
	//Extension of Birk attenuation law proposed by Chou
	// see G.V. O'Rielly et al. Nucl. Instr and Meth A368(1996)745
	//
	double C = 9.59 * 1E-4 * mm * mm / MeV / MeV;
	double response = destep;
	if (birks * destep * stepl * charge != 0.) {
		response = destep / (1. + birks * destep / stepl + C * pow(destep / stepl, 2.));
	}
	return response;
}

map<string, vector<int> > crs_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map<string, vector<int> > MH;

	return MH;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> crs_HitProcess::electronicNoise() {
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
map<int, vector<double> > crs_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map<int, vector<double> > CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double crs_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

