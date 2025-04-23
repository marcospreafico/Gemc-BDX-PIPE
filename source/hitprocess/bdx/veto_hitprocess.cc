// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "veto_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

double att_z_IV = 350 * cm;
double att_z_OV = 350 * cm;

map<string, double> veto_HitProcess::integrateDgt(MHit* aHit, int hitn) {
    map<string, double> dgtz;
    vector<identifier> identity = aHit->GetId();
    
    // channel identifiers
    int sector = identity[0].id; // module number
    int veto_id = identity[1].id; // veto type (now always 4)
    int channel = identity[2].id; // channel
    // 3 front, 4 back, 6 left , 5 right, 1 top, 2 bottom
    int mod = identity[3].id; // module again
    
    // Digitization Parameters
    
    // here all parameters are reported for each fiber
    double optical_coupling[4] = {1, 1, 1, 1};
    double light_yield[4] = {9200, 9200, 9200, 9200};
    double att_length_plastic[4] = {100 * cm, 100 * cm, 100 * cm, 100 * cm};
    double att_length_fiber[4] = {100 * cm, 100 * cm, 100 * cm, 100 * cm};
    double veff = 13 * cm / ns; // this is one for all fibers I guess
    
    // ideally here I could make something as fancy as with crystals where I can specify the SiPM type to make comparisons
    double sensor_surface = pow(0.3*cm, 2);
    double sipm_pde = 0.5;
    
    double length, sside, lside;
    if(channel == 5 || channel == 6){ // left or right
        length = aHit->GetDetector().dimensions[2]*2, sside = aHit->GetDetector().dimensions[0]*2, lside = aHit->GetDetector().dimensions[1]*2;
    }
    if(channel == 1 || channel == 2){ // top or bottom
        length = aHit->GetDetector().dimensions[2]*2, sside = aHit->GetDetector().dimensions[1]*2, lside = aHit->GetDetector().dimensions[0]*2;
    }
    if(channel == 3 || channel == 4){ // front or back
        length = aHit->GetDetector().dimensions[1]*2, sside = aHit->GetDetector().dimensions[2]*2, lside = aHit->GetDetector().dimensions[0]*2;
    }
    
    // readout position
    double x_readout[4] = {0., 0., 0., 0.};
    // here I assume that readout is evenly spaced
    for(int ii = 0; ii < 4; ii++){
        x_readout[ii] = (2.*ii+1.)/8.*lside-lside/2;
    }
    double y_readout = length;
    
    // WARNING - this value has to be updated
    double light_coll[4] = {1., 1., 1., 1.}; // fraction of light collected (for crystal is the ratio of surfaces
    
    // initialize ADC and TDC
    double tdc_conv_crs = 1. / ns;               // TDC conversion factor
    double T_offset_crs = 0 * ns;
    
    double etot[4] = {0., 0., 0., 0.};
    double time[4] = {0., 0., 0., 0.};
    double ADC[4] = {0, 0, 0, 0};
    double TDC[4] = {4096, 4096, 4096, 4096};
    double xave = 0, yave = 0;
    
    // Get info about detector material to eveluate Birks effect
    double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
    
    vector<G4ThreeVector> Lpos = aHit->GetLPos();
    vector<G4double> Edep = aHit->GetEdep();
    vector<G4double> Dx = aHit->GetDx();
    
    // Charge for each step
    vector<int> charge = aHit->GetCharges();
    vector<G4double> times = aHit->GetTime();
    
    unsigned int nsteps = Edep.size();
    
    double Etot = 0;
    double pe_int[4] = {0., 0., 0., 0.};
    double pe[4] = {0., 0., 0., 0.};
    
    for (unsigned int s = 0; s < nsteps; s++) {
        Etot = Etot + Edep[s];
    }
    
    double Etot_B= 0;
    
    
    if (Etot > 0){
        for (unsigned int s = 0; s < nsteps; s++) {
            double Edep_B = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
            // WARNING birks ignored
            
            Etot_B += Edep_B;
            
            double xhit , yhit;
            if(channel == 5 || channel == 6){ // left or right
                xhit = Lpos[s].y(), yhit = Lpos[s].z();
            }
            if(channel == 1 || channel == 2){ // top or bottom
                xhit = Lpos[s].x(), yhit = Lpos[s].z();
            }
            if(channel == 3 || channel == 4){ // front or back
                xhit = Lpos[s].x(), yhit = Lpos[s].y();
            }
            
            xave += xhit/nsteps; yave += yhit/nsteps;
            
            for(int ii = 0; ii < 4; ii++){
                etot[ii] += Edep_B /2 * exp(-(fabs(xhit-x_readout[ii]))/att_length_plastic[ii]) * exp(-(fabs(yhit-y_readout))/att_length_fiber[ii]);
                //if(sector==0) cout << "Attenuation " << Edep[s] << " " << Edep_B << " " << exp(-(fabs(xhit-x_readout[ii]))/att_length_plastic[ii]) << " " << exp(-(fabs(yhit-y_readout))/att_length_fiber[ii]) << " " <<  fabs(xhit-x_readout[ii]) << " " << fabs(yhit-y_readout) << " " << att_length_plastic[ii] << endl;
                // etot = energy * attenuation due to travel in plastic * attenuation along fiber
                time[ii] += (times[s] + sqrt((xhit-x_readout[ii])*(xhit-x_readout[ii])+(yhit-y_readout)*(yhit-y_readout)) / veff) / nsteps;
            }
        }
        
        // old estimate of signal:
        // ADC1 = G4Poisson(pe_sipm[0] * etot_g4 / 2.05);
        
        for(int ii = 0; ii < 4; ii++){
            pe[ii] = int(etot[ii] * light_yield[ii] * sipm_pde * optical_coupling[ii] * light_coll[ii]); // convert to number of phe
            pe[ii] = int(etot[ii]*20); // measurement: 20 phe/MeV
            pe[ii] = G4Poisson(pe[ii]); // add uncertainty
            pe_int[ii] = pe[ii];
            ADC[ii] = pe[ii]+G4RandGauss::shoot(0.,13.); // add noise like in original sim
            if(ADC[ii]<0) ADC[ii] = 0;
            TDC[ii] = time[ii]; // WARNING no time information here
            if(etot[ii] >0.001){
//                cout << "Total energy on fiber " << ii << " is " << etot[ii] << endl;
//                cout << "Corresponding number of phe on said fiber is " << pe[ii] << endl;
            }
        }
    }
    
    
    //if(sector == 0) cout << Etot << " " << Etot_B << " " << etot[0] << " " << etot[1] << " " << etot[2] << " " << etot[3] << endl;
    
    dgtz["hitn"] = hitn;
    dgtz["sector"] = sector;
    dgtz["veto"] = veto_id;
    dgtz["module"] = mod;
    dgtz["channel"] = channel;
    dgtz["adc1"] = ADC[0];        // output in pe
    dgtz["adc2"] = ADC[1];        //deposited energy in keV
    dgtz["adc3"] = ADC[2];        // ignore
    dgtz["adc4"] = ADC[3];        // ignore
    dgtz["adc5"] = 0;        // ignore
    dgtz["adc6"] = 0;        // ignore
    dgtz["adc7"] = xave;        // ignore
    dgtz["adc8"] = yave;        // ignore
    
    dgtz["tdc1"] = TDC[0];        // output in ps
    dgtz["tdc2"] = TDC[1];        // ignore
    dgtz["tdc3"] = TDC[2];        // ignore
    dgtz["tdc4"] = TDC[3];        // ignore
    dgtz["tdc5"] = 4096;        // ignore
    dgtz["tdc6"] = 4096;        // ignore
    dgtz["tdc7"] = 4096;        // ignore
    dgtz["tdc8"] = 4096;        // ignore
    
    return dgtz;
}

vector<identifier> veto_HitProcess::processID(vector<identifier> id, G4Step *step, detector Detector) {
    id[id.size() - 1].id_sharing = 1;
    return id;
}


typedef struct groupedHits {
    double Esum;
    double avgT;
    G4ThreeVector avgPos;
    G4ThreeVector avgLPos;
};

typedef struct aStep {
    bool isMatchedToGroupedHit;
    double E, T;
    G4ThreeVector pos;
    G4ThreeVector lpos;
};

//sort in DESCENDING energy
bool compareSteps(const aStep &a, const aStep &b) {
    return a.E > b.E;
}

bool isMatched(const aStep &a, const aStep &b, double dT = 1 * ns, double dX = 1 * cm) {
    bool ret = false;
    if ((fabs(a.T - b.T) <= dT) && ((a.pos - b.pos).mag() < dX)) {
        ret = true;
    }
    return ret;
}

double veto_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks) {
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

double veto_HitProcess::BirksAttenuation2(double destep, double stepl, int charge, double birks) {
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

map<string, vector<int> > veto_HitProcess::multiDgt(MHit* aHit, int hitn) {
    map<string, vector<int> > MH;
    
    return MH;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> veto_HitProcess::electronicNoise() {
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
map<int, vector<double> > veto_HitProcess::chargeTime(MHit* aHit, int hitn) {
    map<int, vector<double> > CT;
    
    return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double veto_HitProcess::voltage(double charge, double time, double forTime) {
    return 0.0;
}

