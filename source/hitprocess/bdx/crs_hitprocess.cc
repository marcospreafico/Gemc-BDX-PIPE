 // G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "crs_hitprocess.h"

// Root utils
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> crs_HitProcess::integrateDgt(MHit* aHit, int hitn) {
  
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

    int sector = identity[0].id;    //module
	int xch = identity[1].id;       //coordinates in the module
	int ych = identity[2].id;
    int zch = identity[3].id;
    int SiPM = identity[4].id;
    string crs_material = aHit->GetDetector().GetLogical()->GetMaterial()->GetName();

    // Parameter for Crs
                                            // Babar Crs                        PANDA crystals                      L3 crystals
    map<string, double> light_yield_db = {{"CsI_Tl", 50000 * (1. / MeV)}, {"G4_PbWO4", 310 * (1. / MeV)}, {"G4_BGO", 8000 * (1. / MeV)}};
    //map<string, double> light_yield_db = {{"CsI_Tl", 700 * (1. / MeV)}, {"G4_PbWO4", 310 * (1. / MeV)}, {"G4_BGO", 8000 * (1. / MeV)}};// TMP
    map<string, double> attenuation_lenght_db = {{"CsI_Tl", 30 * cm},     {"G4_PbWO4", 60000 * cm},       {"G4_BGO", 25 * cm}};
    map<string, double> optical_coupling_db = {{"CsI_Tl", 0.6866},        {"G4_PbWO4", 0.9},              {"G4_BGO", 0.7}};
    
    double optical_coupling = optical_coupling_db[crs_material];
    double light_yield_crs = light_yield_db[crs_material];
    light_yield_crs = 600*(1./MeV);
	double att_length_crs = attenuation_lenght_db[crs_material]; // compatible with NO ATT Lenght as measured for cosmic muons
    double veff_crs = 30 / 1.8 * cm / ns;                     // light velocity in crystal
  
    //parameters of SiPM
    int sipm_size = SiPM/100;
    int sipm_cell = SiPM - sipm_size * 100;
    
    double sensor_surface_crs = pow(0.1 * sipm_size * cm, 2);
    
    map<int, double> sipm_pde = {{25, 0.22}, {50, 0.35}, {75, 0.5}}; // 25 and 50 from old code; 75 from hamamatsu datasheet
    double sensor_qe_crs = sipm_pde[sipm_cell]; // 24% sipm 100um 22% sipm 25um 35% sipm 50um
    
    // set parameters for specific crystals
    double readout_surface_crs = 0.;
    double length_crs = 0.;
    if(aHit->GetDetector().type == "Box"){
        double sside_crs = 2 * aHit->GetDetector().dimensions[0];
        double lside_crs = 2 * aHit->GetDetector().dimensions[1];
        readout_surface_crs = sside_crs * lside_crs * mm * mm;
        length_crs = aHit->GetDetector().dimensions[2];
    }else if(aHit->GetDetector().type == "Trd"){
        double sside_crs = 2 * aHit->GetDetector().dimensions[0];
        double lside_crs = 2 * aHit->GetDetector().dimensions[2];
        readout_surface_crs = sside_crs * lside_crs * mm * mm;
        length_crs = aHit->GetDetector().dimensions[4];
    }
	
    
    // matching readout surface to crystal parameters
    
	double light_coll_crs = sensor_surface_crs / readout_surface_crs;
	if (light_coll_crs > 1) light_coll_crs = 1.;
	
	double etotL_crs = 0; //L= Large side redout
	double etotR_crs = 0; //R= short side redout
    
	double timeL_crs = 0;
	double timeR_crs = 0;
	
	double tdc_conv_crs = 1. / ns;               // TDC conversion factor
	double T_offset_crs = 0 * ns;
    
	double ADCL_crs = 0;
	double ADCR_crs = 0;
	double TDCL_crs = 4096;
	double TDCR_crs = 4096;
	
    double TDCB = 4096;

	// Get info about detector material to eveluate Birks effect
	double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
    
    // forcing Birks for CsI(Tl) from arxiv.org/pdf/0911.3041 and checked with alpha
	// birks_constant=3.2e-3 best at  8e-3
    if(crs_material == "CsI_Tl"){
        birks_constant = 3.2e-3; // g / MeV / cm^2
    }else if(crs_material == "G4_PbWO4"){
        birks_constant = 10.5e-3;
    }else if(crs_material == "G4_BGO"){
        // Here the paper does not provide a value for BGO
    }

	double time_min_crs[4] = { 0, 0, 0, 0 };

	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> Dx = aHit->GetDx();

	// Charge for each step
	vector<int> charge = aHit->GetCharges();
	vector<G4double> times = aHit->GetTime();

	unsigned int nsteps = Edep.size();
    
	double Etot = 0;

	double* test;
	double tim;
    
	double peR_int_crs = 0;
	double peR_crs = 0.;
    
	double peL_int_crs = 0;
	double peL_crs = 0.;
    
    int Nsamp_int = 500; // 2 us (BDX proto) time window
    // here you can change integration time depending on the material
    if(crs_material == "CsI_Tl"){
        Nsamp_int = 500;
    }else if(crs_material == "G4_PbWO4"){
        Nsamp_int = 500;
    }
	double sigmaTR_crs = 0.; //??
    
	for (unsigned int s = 0; s < nsteps; s++) {
		Etot = Etot + Edep[s];
	}

	double Etot_B= 0;
    
	if (Etot > 0) {
		for (unsigned int s = 0; s < nsteps; s++) {   //Reference vie for cal matrix:
													  //cristals with short size pointing downstream
													  // sipm attached to the large side (upstream)
													  // left: smoll size, right: large size
													  // Use only dRight
													  // for rotated (old) crystal we keep the same convention:
													  // readout = small size (use dLeft)
            
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
        //cout << peR_crs << endl;
		peR_crs = G4Poisson(peR_crs);
        
        //WARNING
        test = WaveForm(peR_crs, &tim, crs_material);
        
        TGraph* WF_new = new TGraph();
        for(int s = 0; s < Nsamp_int; s++){
            WF_new->SetPoint(s, s, test[s]);
        }
        
        TFile* fout = new TFile("testWF.root", "RECREATE");
        fout->cd();
        WF_new->Write(); fout->Write(); fout->Close();
        
        double peR_int_crs_old = 0;
		for (unsigned int s = 0; s < Nsamp_int; s++) {
			peR_int_crs += test[s];
		}
        
        //      Left readout
        peL_crs = int(etotL_crs * light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs);
        peL_crs = G4Poisson(peL_crs);
        test = WaveForm(peL_crs, &tim, crs_material);
        
        for (unsigned int s = 0; s < Nsamp_int; s++) {
            peL_int_crs = peL_int_crs + test[s];
        }
        
        peL_int_crs = peL_crs;
        peR_int_crs = peR_crs;
        
        // Save variables
        
        // Correct readout to be in energy;
        // This variable keeps track of the fraction of energy measured in the integration window and corrects for it
//        double  ts =  0.680, fs =  0.64, tl =  3.34, fl = 0.36; // fraction of long /short time; value of long/short time
//        if(crs_material == "G4_PbWO4"){
//            ts = 0.00680; fs= 0.64; tl = 0.0334; fl =  0.36;
//        } else if(crs_material == "CsI_Tl"){
//            ts =  0.680; fs =  0.64; tl =  3.34; fl = 0.36;
//        }
//        double digiframe =  Nsamp_int * 4. / 1000.;
//        double sigfrac = 1 - (fs* exp(-digiframe / ts) + fl * exp(-digiframe / tl)); // fraction of signal contained in a digiframe digitalization window
        
        // energy measured = numbrer of phe / (light yield * attenuation * light was 2x before splitting * fraction measured in a fixed time window
        ADCR_crs = (peR_int_crs)/(light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs * 0.5); // * sigfrac); // in MeV
        ADCL_crs = (peL_int_crs)/(light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs * 0.5); // * sigfrac);
        
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
	dgtz["adcf"] = ADCL_crs * 1000000;
	dgtz["tdcb"] = TDCB * 1000.;	  //original time in ps
	dgtz["tdcf"] = Etot*1000000;
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

double* crs_HitProcess::WaveForm(double npe, double* time, string crs_material){
    double c = exp(-2);
    
    int Nsamp_WF = 2500; // number of samples to generate WF; 2500 = 10 us
    
    static double* WFsample = new double[Nsamp_WF]; // object to save the WF
    for(unsigned int s = 0; s < Nsamp_WF; s++){ WFsample[s] = 0; }
    
    double smp_t = 4./1000.; // assuming fADC sampling at 250 MHz 1 sample = 4 ns
    
    array<double, 6> p = {0., 0., 0., 0., 0., 0.}; // crs scintillation time parameters
    if(crs_material == "G4_PbWO4"){
        p = { 0., 0.00680, 0.64, 0.0334, 0.36, 0. }; // PbWO4
    }
    else if(crs_material == "CsI_Tl"){
        p = { 0., 0.680, 0.64, 3.34, 0.36, 0. }; // BaBar CsI
    }
    else{
        cout << "ERROR: no waveform parametrization for this material ( "+crs_material+" )" << endl;
        return 0;
    }
    TF1* tdistrib = new TF1("tdistrib", "([2]/[1]*exp(-x/[1])+[4]/[3]*exp(-x/[3]))/([2]/[1]+[4]/[3])", 0, Nsamp_WF*smp_t);
    for(int ii = 0; ii < 6; ii ++){ tdistrib->SetParameter(ii, p[ii]); }
    
    
    //Definition of phe shape
    double tau = 15.; // ampli response time constant (in ns)
    double t0 = 0.01; // t0 starting time (in ns)
    double area = (tau / c / 2.);
    double A = 1. / area; // amplitude at mnax (55.41 to have it normalized to integral=1, otherwise the max is at 1)
    
    double t_spread = 1. * 0.000; // pream time spread in us
    double A_spread = 1. * 0.4 * A; // pream amp spread (in fraction of 1pe amplitude = A)
    if(crs_material == "CsI_Tl"){
         t_spread = 1. * 0.020; // pream time spread in us
         A_spread = 1. * 0.05 * A; // pream amp spread (in fraction of 1pe amplitude = A)
    }
    
    static double AmpWF[80]; // phe WF
    for (unsigned int s = 0; s < 80; s++) {
        double t = 1000. * s * smp_t;
        double func = (t - t0) * (t - t0) * exp(-(t - t0) / tau) * A / (4 * tau * tau * c) * 0.5 * (abs(t - t0) / (t - t0) + 1);
        AmpWF[s] = smp_t * 1000. * func;
    }
    
    static double frac = 1 - ((p[2] * exp(-smp_t * Nsamp_WF / p[1]) + p[4] * exp(-smp_t * Nsamp_WF / p[3])));    // fraction of pe in Nsamp_WF
    
    
    // generate waveform sample
    double t; int it; // time variable
    int mNpe = G4Poisson(frac * npe); // number of phe in Nsamp_WF
    for(unsigned int s = 1; s <= mNpe; s++){
        t = tdistrib->GetRandom();
        t = G4RandGauss::shoot(t, t_spread);
        if(t < 0.) t = 0.;
        if(t > smp_t * Nsamp_WF) t = smp_t * Nsamp_WF;
        it = t / smp_t;
        for(unsigned int s = 0; s < 80; s++){ // sum the phe WF at the phe time
            double func = G4RandGauss::shoot(AmpWF[s], A_spread);
            if((it + s) > Nsamp_WF) break;
            WFsample[it + s] += func;
        }
    }
    
    // mimicking a CF discriminator at 1/3 of the max signal
    *time = 0.;
    double time_max = -100;
    int s = 0, s_time_max = 0;
    while(time_max < WFsample[s]){// search for max WF time
        time_max = 1/2. * (WFsample[s + 1] + WFsample[s]);
        s_time_max = s;
        *time = 1000. * smp_t * s_time_max / 3.;
        s++;
    }
    
    return WFsample;
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

