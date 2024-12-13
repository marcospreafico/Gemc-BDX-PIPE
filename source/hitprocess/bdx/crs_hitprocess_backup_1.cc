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
    string crs_material = aHit->GetDetector().GetLogical()->GetMaterial()->GetName();

    // Parameter for BaBar Crystal
    double optical_coupling = 0.95;
	double light_yield_crs = 50000 * (1. / MeV);
	double att_length_crs = 60 * cm; // compatible with NO ATT Lenght as measured for cosmic muons
	double sensor_surface_crs = pow(0.6 * cm, 2);
	double redout_surface_crs = pow(5.1 * cm, 2);
	double sensor_qe_crs = 0.22; // 24% sipm 100um 22% sipm 25um 35% sipm 50um
    double veff_crs = 30 / 1.8 * cm / ns;                     // light velocity in crystal
    
    // set parameters for specific crystals
    if(crs_material == "CsI_Tl"){
        sensor_surface_crs = pow(0.6 * cm, 2);
        sensor_qe_crs = 0.22; // consider only 25um sipm
        optical_coupling = 0.6866;
        att_length_crs = 60 * cm;
        light_yield_crs = 50000 * (1. / MeV);
    }else if(crs_material == "G4_PbWO4"){
        sensor_surface_crs = pow(0.6 * cm, 2);
        sensor_qe_crs = 0.22; // consider only 25um sipm
        optical_coupling = 0.9;
        att_length_crs = 60000 * cm; // compatible with NO ATT Lenght as measured for cosmic muons
        light_yield_crs = 310 * (1. / MeV); //Panda Crystals LY
    }else if(crs_material == "G4_BGO"){
        sensor_surface_crs = pow(0.3 * cm, 2);
        sensor_qe_crs = 0.5; // source: Tommaso
        optical_coupling = 0.7;
        att_length_crs = 25 * cm;
        light_yield_crs = 8000 * (1. / MeV);
    }
    
	double length_crs = 2 * aHit->GetDetector().dimensions[4]; //(in mm)
	double sside_crs = 2 * aHit->GetDetector().dimensions[0];
	double lside_crs = 2 * aHit->GetDetector().dimensions[2];
    
    // matching readout surface to crystal parameters
	redout_surface_crs = sside_crs * lside_crs * mm * mm;

	double light_coll_crs = sensor_surface_crs / redout_surface_crs;
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
		peR_crs = G4Poisson(peR_crs);
        
        //WARNING
        test = WaveForm(peR_crs, &tim, crs_material);
        double * test_old = WaveFormPbwo(peR_crs, &tim);
        
        //DEBUG: plot WF
        TGraph* WF_old = new TGraph();
        TGraph* WF_new = new TGraph();
        for(int s = 0; s < Nsamp_int; s++){
            WF_old->SetPoint(s, s, test_old[s]);
            WF_new->SetPoint(s, s, test[s]);
        }
        TFile* fout = new TFile("testWF.root", "RECREATE");
        fout->cd();
        WF_old->Write(); WF_new->Write(); fout->Write(); fout->Close();
        
        
//        //if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "CsI_Tl"){
//            test = WaveForm(peR_crs, &tim);
//        //}else
//        if(aHit->GetDetector().GetLogical()->GetMaterial()->GetName() == "G4_PbWO4"){
//            test = WaveFormPbwo(peR_crs, &tim);
//        }
        double peR_int_crs_old = 0;
		for (unsigned int s = 0; s < Nsamp_int; s++) {
			peR_int_crs += test[s];
            peR_int_crs_old += test_old[s];
		}
        
        cout << peR_int_crs / peR_int_crs_old << " " << peR_int_crs << " " << peR_int_crs_old << endl;
        
        fold << peR_int_crs_old << endl; fnew << peR_int_crs << endl;
        
        
        //      Left readout
        peL_crs = int(etotL_crs * light_yield_crs * sensor_qe_crs * optical_coupling * light_coll_crs);
        peL_crs = G4Poisson(peL_crs);
        if(crs_material == "G4_PbWO4"){
            test = WaveFormPbwo(peL_crs, &tim);
        } else if(crs_material == "CsI_Tl"){
            test = WaveForm(peL_crs, &tim, crs_material);
        }
        
        for (unsigned int s = 0; s < Nsamp_int; s++) {
            peL_int_crs = peL_int_crs + test[s];
        }
        
        // Save variables
        
        // Correct readout to be in energydouble sigfrac = 0;
        double  ts =  0.680, fs =  0.64, tl =  3.34, fl = 0.36; // fraction of long /short time; value of long/short time
        if(crs_material == "G4_PbWO4"){
            ts = 0.00680; fs= 0.64; tl = 0.0334; fl =  0.36;
        } else if(crs_material == "CsI_Tl"){
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

//double crs_HitProcess::WaveForm(double npe, double time)
double* crs_HitProcess::WaveForm_old(double npe, double* time) {
	double c = exp(-2.);
	int it;
    
    // initialization of waveform
	int Nch_digi = 2500; //Number of samples in waveform - 2500 = 10 us
    static double* WFsample = new double[Nch_digi+200]; // Needs to be >  Nch_digi+size of the response to the single pe
    for (unsigned int s = 0; s < Nch_digi+200; s++) { WFsample[s] = 0; }
	double smp_t = 4. / 1000.; // Assuming fADC sampling at 250 MHz 1sample every 4ns

    // parameters of time distribution
	double p[6] = { 0., 0.680, 0.64, 3.34, 0.36, 0. }; // Babar CsI paprameters: 0, fast component(us), % fast, slow comp(in us), % slow, 0
    TF1* tdistrib = new TF1("tdistrib", "([2]/[1]*exp(-x/[1])+[4]/[3]*exp(-x/[3]))/([2]/[1]+[4]/[3])", 0, Nch_digi*smp_t);
    for(int ii = 0; ii < 6; ii ++){ tdistrib->SetParameter(ii, p[ii]); }

    // parameters of the signal waveform
	double tau = 15.; // ampli response time constant (in ns)
	double t0 = 0.01; // t0 starting time (in ns)
	double area = (tau / c / 2.);
	double A = 1. / area; // amplitude at mnax (55.41 to have it normalized to integral=1, otherwise the max is at 1)
    
    // spreads
	double t_spread = 1. * 0.020; // pream time spread in us
	double A_spread = 1. * 0.05 * A; // pream amp spread (in fraction of 1pe amplitude = A)
    
    
	// Building the response to a single pe (preamps response)
    static double AmpWF[80];
    for (unsigned int s = 0; s < 80; s++) {
        double t = 1000. * s * smp_t;
        double func = (t - t0) * (t - t0) * exp(-(t - t0) / tau) * A / (4 * tau * tau * c) * 0.5 * (abs(t - t0) / (t - t0) + 1);
        AmpWF[s] = smp_t * 1000. * func;
    }
    static double frac = 1 - ((p[2] * exp(-smp_t * Nch_digi / p[1]) + p[4] * exp(-smp_t * Nch_digi / p[3])));// fraction of pe in Nch_digi
  
    npe = npe*frac;
   
    /*for (unsigned int s = 1; s <= npe; s++) {
        double t = tdistrib->GetRandom(0, Nch_digi * smp_t); // time in usec
        // spreading time and amplitude of the ampli signal
        t = G4RandGauss::shoot(t, t_spread);
        if (t < 0.) t = 0.;
        it = t / smp_t;

        for (unsigned int s = 0; s < 80; s++) {
            t = 1000. * s * smp_t;
            double func = AmpWF[s];
            func = G4RandGauss::shoot(func, A_spread);
            if ((s + it) < Nch_digi) WFsample[s + it] = WFsample[s + it] + func;
        }
    }*/
    
    double y, rr, WF;
    for (unsigned int s = 1; s <= npe; s++) {
        y = 1.;
        WF = 0.;
        double t;
        while (y > WF) {
            rr = (rand() % 1000000 + 1) / 1000000.; // rnd number between 0-1
            t = Nch_digi * smp_t * rr; // extracting over 5000 samples range (5000x4ns=20us)
            //WF= 1./5.15*((1-exp(p[0]+p[1]*t))*exp(p[2]+p[3]*t)+exp(p[4]+p[5]*t));
            WF = (p[2] / p[1] * exp(-t / p[1]) + p[4] / p[3] * exp(-t / p[3])) / (p[2] / p[1] + p[4] / p[3]); //pulire facendo estrazione da questa funzione
            rr = (rand() % 1000000 + 1) / 1000000.; // rnd number between 0-1
            y = rr;
            //  cout << "WF " << WF   << " rnd " << y  << endl;
        }
        // spreading time and amplitude of the ampli signal
        t = G4RandGauss::shoot(t, t_spread);
        if (t < 0.) t = 0.;
        it = t / smp_t;

        for (unsigned int s = 0; s < 80; s++) {
            t = 1000. * s * smp_t;
            double func = AmpWF[s];
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

	return WFsample;

}

double* crs_HitProcess::WaveFormPbwo(double npe, double* time_pbwo) {
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

