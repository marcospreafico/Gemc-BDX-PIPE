
// G4 headers
#include "JPOS_HCAL_hitprocess.h"

#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> JPOS_HCAL_HitProcess::integrateDgt(MHit* aHit, int hitn) {
  
  map<string, double> dgtz;

  vector<identifier> identity = aHit->GetId(); 
  int sector = identity[0].id;
  int layer = identity[1].id;
  int channel = identity[2].id;
  //cout<<"sector  "<<sector<<endl;
  //cout<<"layer  "<<layer<<endl;
  //cout<<"channel  "<<channel<<endl;
		       

  double d0 = aHit->GetDetector().dimensions[0]; //semidimensioni sbarretta in mm risentono delle rotazioni 
  double d1 = aHit->GetDetector().dimensions[1];
  double d2 = aHit->GetDetector().dimensions[2];
  //cout<<"dimensions "<<d0<<" "<<d1<<" "<<d2<<endl;


  // Get info about detector material to eveluate Birks effect
  double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
  //cout << "Birks constant is: " << birks_constant << endl;
  //cout << aHit->GetDetector().GetLogical()->GetMaterial()->GetName() << endl;


  vector<G4double> Edep = aHit->GetEdep();//indicizzato sul numero di step  nel rilevatore
  vector<G4double> times = aHit->GetTime();
  vector<G4ThreeVector> Lpos = aHit->GetLPos();
  vector<G4double> Dx = aHit->GetDx();
  vector<int> charge = aHit->GetCharges();

  int nsteps = Edep.size();
  //cout<<"hitn "<< hitn<<endl;
  //cout<<"n steps "<< nsteps<<endl;

  double t0 = times[0];
  
 
  double Etot = 0;
  double EtotB = 0;

  double Mx = 0;
  double My = 0;
  double Mz = 0;
	
  for (int s = 0; s < nsteps ; s++) {

    Etot += Edep[s];

    double Edep_B = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
    EtotB += Edep_B;
    //cout<<Lpos[s].getX()<<endl;                 //posizione step relativa alla sbarretta
    Mx += Lpos[s].getX();
    //cout<<Lpos[s].getY()<<endl;
    My += Lpos[s].getY();
    //cout<<Lpos[s].getZ()<<endl;
    Mz += Lpos[s].getZ();
    
  }

  Mx = Mx/nsteps;
  My = My/nsteps; 
  Mz = Mz/nsteps;

  double hit_pos=4004;  
  if((d0>d1)&&(d0>d2)){
    hit_pos=Mx;
    //cout<<"x"<<endl;
  }
  else if((d1>d0)&&(d1>d2)){
    hit_pos=My;
    //cout<<"y"<<endl;
  }
  else if((d2>d0)&&(d2>d1)){
    hit_pos=Mz;
    //cout<<"z"<<endl;
  }

  //cout<<"Etot "<<Etot<<" EtotB "<<EtotB<<endl;
  //cout<<"hit_pos "<<hit_pos<<endl;
  //cout<<endl;
  	      
  
  
  dgtz["hitn"] = hitn;
  dgtz["sector"] = sector;
  dgtz["layer"] = layer;
  dgtz["channel"] = channel;
  dgtz["t0"] = t0;           // ns
  dgtz["Etot"] = Etot;       // MeV
  dgtz["EtotB"] = EtotB;     // Deposited energy corrected by Birks' Law MeV
  dgtz["hit_pos"] = hit_pos; // mm

  return dgtz;
}




double JPOS_HCAL_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks) {
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

double JPOS_HCAL_HitProcess::BirksAttenuation2(double destep, double stepl, int charge, double birks) {
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




vector<identifier> JPOS_HCAL_HitProcess::processID(vector<identifier> id, G4Step *step, detector Detector) {
	id[id.size() - 1].id_sharing = 1;
	return id;
}


map<string, vector<int> > JPOS_HCAL_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map<string, vector<int> > MH;

	return MH;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> JPOS_HCAL_HitProcess::electronicNoise() {
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
map<int, vector<double> > JPOS_HCAL_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map<int, vector<double> > CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double JPOS_HCAL_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

