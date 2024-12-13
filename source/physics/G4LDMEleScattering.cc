/*
 * G4LDMEleScattering.h
 *
 *  Created on: March 18th, 2024
 *      Author: mspreafico
 */

#include "G4LDMEleScattering.h"

#include <iostream>
#include <map>
#include <stdlib.h>
#include <sys/stat.h>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4EmProcessSubType.hh"
#include "G4NistManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4PrimaryParticle.hh"

#include "G4DarkMatter.h"
#include "G4SystemOfUnits.hh"

//root headers
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
#include "TObject.h"
#include "TList.h"
#include "TKey.h"
#include "TRandom.h"

//Keeping global variables shared possibly by different istances of the same file means this is not thread-safe.
//Will use gemc in single-thread mode.

static bool messageOut = true;

//This holds the total cross section, as a function of the impinging DM total energy
static TGraph* xsecData;

G4LDMEleScattering::G4LDMEleScattering(const G4String& processName, G4ProcessType type) :
    G4VDiscreteProcess(processName, type), MDM(G4DarkMatter::DarkMatter()->GetPDGMass()), LowestEnergyLimit(1000 * eV), //the lower limit of the minimum energy
    HighestEnergyLimit(1e21 * eV), // ok to 1e21eV=1e12GeV, then LPM suppression
    CrossSecFactor(1.), fData(0), gXsecData(0), hXsecDistr(0) { // !! Here I set the CrossSecFactor to be eps/eps0 with eps0 = 1e-3
        SetProcessSubType(G4LDMEleScatteringProcessID);
        MeanFreePath = DBL_MAX;

	/*The data is loaded "on the fly", in the GetMeanFreePath - for the given ELEMENT.*/
	char* pPath;
        G4double LDMSpin = 0.5;
        if(LDMSpin == 0.5) pPath = getenv("FERMION_LDM");
        else if(LDMSpin == 0.) pPath = getenv("SCALAR_LDM");
        else{
            std::cout << "ERROR: particle spin not compelling with models considered (spin = 0, 1/2). EXITING" << std::endl;
            exit(1);
        }
	if (pPath == NULL) {
		std::cout << std::endl << std::endl;
		std::cout << "Environmental variable for database FERMION_LDM or SCALAR_LDM is not set." << std::endl;
		std::cout << "This should be set to the folder with the data for the LDM interaction" << std::endl;
		std::cout << "The folder should contain cross sections in file named data_DarkPhotonMAssInGeV.root" << std::endl;
		std::cout << "The DarkPhotonMAssInGeV should have a precision of 3 (i.e. 0.010 is ok, 0.01 is not for 10 MeV)" << std::endl;
		exit(1);
	} else {
		dataPathFolder = G4String(pPath);
		std::cout << std::endl << std::endl;
		if (messageOut) {
			std::cout << "G4ProcessSkeleton creator, the environmental variable " << pPath <<  " is set to: " << dataPathFolder << std::endl;
			messageOut = false;
		}
		dataPathFolder += "/LDM"; // in the original scalar process here I was taking the /muon folder --> mother particles 2
		struct stat sb;
		if (!(stat(dataPathFolder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) {
			std::cout << std::endl << std::endl;
			std::cout << "G4ParticleSkeleton, LDM folder is missing within CROSS_SECTION_PATH" << std::endl;
			exit(1);
		}

	}
        
        gRandom->SetSeed(0);
}

G4LDMEleScattering::~G4LDMEleScattering() {
	// TODO Auto-generated destructor stub
}

G4bool G4LDMEleScattering::IsApplicable(const G4ParticleDefinition& particle) {
	return (&particle == G4DarkMatter::DarkMatter());
}

G4double G4LDMEleScattering::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*) {

// returns the particle mean free path in GEANT4 internal units
// (MeanFreePath is a private member of the class)

	aTrack.GetParticleDefinition()->GetPDGMass();

	const G4DynamicParticle* aDynamicParticle = aTrack.GetDynamicParticle();
	G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();
	//G4double Mass = aTrack.GetParticleDefinition()->GetPDGMass();
	G4Material* aMaterial = aTrack.GetMaterial();

	if (KineticEnergy <= LowestEnergyLimit) MeanFreePath = DBL_MAX;
	else
		MeanFreePath = ComputeMeanFreePath(KineticEnergy, aMaterial);

	return MeanFreePath;
}

G4double G4LDMEleScattering::ComputeMeanFreePath(G4double KineticEnergy, G4Material* aMaterial) {

// computes and returns the photon mean free path in GEANT4 internal units

	const G4ElementVector* theElementVector = aMaterial->GetElementVector();
	const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

	G4double nSIGMA = 0;

    for (size_t i = 0; i < aMaterial->GetNumberOfElements(); ++i) {
        //G4double AtomicZ = (*theElementVector)[i]->GetZ();
        //G4double AtomicA = (*theElementVector)[i]->GetA() / (g / mole);
        G4String AtomicSymbol = (*theElementVector)[i]->GetSymbol();
        //nSIGMA += NbOfAtomsPerVolume[i] * ComputeCrossSectionPerAtom(KineticEnergy, AtomicSymbol);
        nSIGMA += NbOfAtomsPerVolume[i] * ComputeCrossSectionPerAtom(KineticEnergy, (*theElementVector)[i]);
    }
	return nSIGMA > DBL_MIN ? 1. / nSIGMA : DBL_MAX;
}

G4double G4LDMEleScattering::GetCrossSectionPerAtom(const G4DynamicParticle* aDynamicParticle, G4Element* anElement)

// gives the total cross section per atom in GEANT4 internal units
		{
	G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();
	G4String AtomicSymbol = anElement->GetSymbol();
    G4double AtomicZ = anElement->GetZ();
	//G4double crossSection = AtomicZ * ComputeCrossSectionPerAtom(KineticEnergy, AtomicSymbol);
            G4double crossSection = AtomicZ * ComputeCrossSectionPerAtom(KineticEnergy, anElement);

            //            // !! Here I am looking for DM scattering off electrons so I need to scale for the number of electrons in the nucleus
            //std::cout << " Cross section for " << AtomicSymbol << " is " << crossSection << std::endl;
	return crossSection;
}

//G4double G4LDMEleScattering::ComputeCrossSectionPerAtom(G4double Ekin, G4String ElementSymbol) {
G4double G4LDMEleScattering::ComputeCrossSectionPerAtom(G4double Ekin, G4Element* anElement) {

	// Calculates the microscopic cross section in GEANT4 internal units.
	// Note that this function also reads data "on the fly", storing data in the map

	if (Ekin <= LowestEnergyLimit) return 0; // below threshold return 0

    // if the xsecValData does not exist, I read the files to find the total cross section graph
    if (xsecData == NULL) {
        fData = this->GetDataFile();
        if (fData->IsZombie()) {
            std::cout << "File does not exists! Return 0 xsection" << std::endl;
            delete fData;
            return 0.0;
        }
        
        gXsecData = NULL;
        gXsecData = (TGraph*) fData->Get("gData");
        if (gXsecData == NULL) {
            std::cout << "Data does not exists! Return 0 xsection" << std::endl;
            delete fData;
            return 0.0;
        }
        /*I clone the original data graph, then delete the file*/
        xsecData = (TGraph*) gXsecData->Clone("gData");
        delete fData;
    }
    
	G4double CrossSection = 0.0;

	CrossSection = gXsecData->Eval(Ekin/GeV); //Need to evaluate in GeV!!!
    //std::cout << Ekin/GeV << std::endl;
    G4String AtomicSymbol = anElement->GetSymbol();
    G4double AtomicZ = anElement->GetZ();
    CrossSection =  AtomicZ * CrossSection * picobarn; //data is in pbarn
        
    //std::cout << "Cs from graph: " << picobarn << std::endl;
    //std::cout << "Element: " << AtomicSymbol << " with Z = " << AtomicZ << std::endl;
    
    //std::cout << " Cross section for " << AtomicSymbol << " is " << CrossSection/100. << std::endl;
    //CrossSecFactor = eps/1e-3; // !! Here I scale the cross section factor to the value of epsilon I want to use
    //std :: cout << CrossSecFactor << std::endl;
	CrossSection *= CrossSecFactor; // increase the CrossSection by  (by default 1)
	return CrossSection;
}

G4VParticleChange* G4LDMEleScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
//
// interaction between LDM and detecotr
//
  //std::cout << "Scalar production" << std::endl;

	G4double deltaMin = DBL_MAX;
	G4double delta;
	G4int idx;

	G4double randomE, randomZ, randomPhi;

	aParticleChange.Initialize(aTrack);
	G4Material* aMaterial = aTrack.GetMaterial();

	// current LDM kinetic energy and direction, return if energy too low
	const G4DynamicParticle *aDynamicParticle = aTrack.GetDynamicParticle();
	G4double kineticE = aDynamicParticle->GetKineticEnergy();
	if (kineticE <= LowestEnergyLimit) {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}
    
	G4ParticleMomentum MomentumDirection = aDynamicParticle->GetMomentumDirection();

	// select randomly one element constituting the material
	const G4Element* anElement = SelectRandomAtom(aDynamicParticle, aMaterial);
	//try again
	if (anElement == NULL){
		 anElement = SelectRandomAtom(aDynamicParticle, aMaterial);
	}
	//if still bad, do nothing
	if (anElement == NULL){
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	//Get the graph with the data
	gXsecData = xsecData;
    
    G4double EDM = kineticE + MDM; // !! Here I need to switch to the total energy instead of the kinetic energy

    //Compute the index for this kinetic energy E -> the index of the closest point. I assume nothing about the points (not equally spaced, nor ordered)
	for (int ip = 0; ip < gXsecData->GetN(); ip++) {
    	delta = fabs(kineticE / GeV - gXsecData->GetX()[ip]);
		if (delta < deltaMin) {
			deltaMin = delta;
			idx = ip;
		}
	}

	fData = this->GetDataFile();
	if (fData->IsZombie()) {
		std::cout << "File does not exists! Return 0 xsection" << std::endl;
		delete fData;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}
	hXsecDistr = (TH2D*) fData->Get(Form("h2_%i", idx));

	//Following error may be due to binning+threshold effects
	if (hXsecDistr == NULL) {
		std::cout << "Error! Cannot find histo. Doing nothing. Ekin: " << kineticE << " idx: " << idx << std::endl;
		delete fData;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	/*If we are here, it means we have the histogram. Get random number from it: fraction of total energy of impinging DM particle inherited by electron
     and cosine of scattering angle (Y)*/
	hXsecDistr->GetRandom2(randomE, randomZ);
	int repeat_extraction = 0; 
    
    G4double me = G4Electron::Electron()->GetPDGMass();
    
    // I repeat the extraction if the electron inherits too much energy and I would end up with the DM particle having E < M_DM
	while(((EDM+me-randomE*EDM) < MDM || randomE*EDM < me) && repeat_extraction<10){
	  hXsecDistr->GetRandom2(randomE, randomZ);
	  std::cout << "Final state DM with negative kinetic energy" << std::endl;
	  repeat_extraction ++;
	}
	randomE = randomE * EDM; //always need to specify the units in G4.
	if(repeat_extraction == 9) {
	  std::cout << "End extraction, could not find correct combination" << std::endl; 
	}

	//Just in case..
	if (randomZ > 1) randomZ = 1;
	if (randomZ < -1) randomZ = -1;
	randomPhi = G4UniformRand()*twopi;
	G4ThreeVector ElectronDirection(sqrt(1 - randomZ * randomZ) * sin(randomPhi), sqrt(1 - randomZ * randomZ) * cos(randomPhi), randomZ);
	// rotate to actual impinging muon direction
	ElectronDirection.rotateUz(MomentumDirection);

	/*Define the scattered muon momentum as the difference: impinging - scalar*/
	G4ThreeVector ImpingingMomentum = aDynamicParticle->GetMomentum();
	G4ThreeVector ElectronMomentum = ElectronDirection * sqrt(randomE * randomE - MDM * MDM);
	G4ThreeVector DMMomentum = ImpingingMomentum - ElectronMomentum;

    std::ofstream ofile("test.txt", std::ios_base::app);
    std::cout << "DM interaction" << std::endl;
    std::cout << "Impinging total energy: " << EDM << "; final electron energy " << randomE << "; scattering cosine " << randomZ << std::endl;
    ofile << EDM << " " << randomE << " " << randomZ << std::endl;
    
	aParticleChange.SetNumberOfSecondaries(2);
	G4DynamicParticle* aParticle1 = new G4DynamicParticle(G4Electron::Electron(), ElectronDirection, randomE);
	aParticleChange.AddSecondary(aParticle1);
    
	// Change the incident particle
	/*aParticleChange.ProposeMomentumDirection(DMMomentum.unit());
	G4double finalKE = EDM+me-randomE;
	if(finalKE<0){
	  std::cout << "WARNING: negative final kinetic energy - " << finalKE << std:: endl;
	  finalKE = 0; 
	}
	aParticleChange.ProposeEnergy(EDM + me - randomE - MDM); // E' = E0 - Ee - MDM (here is a kinetic energy)
	std::cout << kineticE << std::endl;
     */
	aParticleChange.ProposeTrackStatus(fStopAndKill);

	//  Reset NbOfInteractionLengthLeft and return aParticleChange
	delete fData;


	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4Element* G4LDMEleScattering::SelectRandomAtom(const G4DynamicParticle* aDynamicParticle, G4Material* aMaterial) {
	// select randomly 1 element within the material, invoked by PostStepDoIt

	const G4int NumberOfElements = aMaterial->GetNumberOfElements();
	const G4ElementVector* theElementVector = aMaterial->GetElementVector();
	if (NumberOfElements == 1) return (*theElementVector)[0];
	
	const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

	G4double PartialSumSigma = 0.;
	G4double rval = G4UniformRand()/MeanFreePath;

	for (G4int i = 0; i < NumberOfElements; ++i) {
		PartialSumSigma += NbOfAtomsPerVolume[i] * GetCrossSectionPerAtom(aDynamicParticle, (*theElementVector)[i]);
		//std::cout<<"SelectRandomAtom "<<rval<<" "<<PartialSumSigma<<" "<<NbOfAtomsPerVolume[i] * GetCrossSectionPerAtom(aDynamicParticle, (*theElementVector)[i])<<std::endl;
		if (rval <= PartialSumSigma) return ((*theElementVector)[i]);
	}
	//std::cout << " WARNING !!! - The Material '" << aMaterial->GetName() << "' has no elements, NULL pointer returned. rval: " <<rval<< std::endl;
	return NULL;
}

void G4LDMEleScattering::SetCrossSecFactor(G4double fac) {
// Set the factor to artificially increase the cross section

	CrossSecFactor = fac;
	std::cout << "The cross section for ScalarProductin is artificially " << "increased by the CrossSecFactor=" << CrossSecFactor << std::endl;
}

/*This retuns a new TFile for the current element symbol.
 * It is up to the user to delete it!
 */
TFile* G4LDMEleScattering::GetDataFile() {
	G4double MAprimeTMP = MDM*3. / GeV; // !! Here I am fixing by hand the mass ratio but I need to take it from the configuration parameters
	G4String massString;
	G4String filePath;

	std::stringstream stream;
	stream << std::fixed << std::setprecision(3) << MAprimeTMP;

	massString = G4String(stream.str().c_str());
	filePath = dataPathFolder + "/data_" + massString + "GeV.root";

	return new TFile(filePath.c_str());
}

void G4LDMEleScattering::SetEps(G4double fac) {
    // Set the factor to artificially increase the cross section
    
    eps = fac;
    std::cout << "The eps value for LDM scattering is set to " << eps << std::endl;
    double eps_0 = 3.87e-4;
    double bias = (eps/eps_0)*(eps/eps_0);
    
    std::cout << "Since the default value of eps is " << eps_0 << " , the cross section has to be scaled by a factor " << bias << std::endl;
    this->SetCrossSecFactor(CrossSecFactor*bias); // here I add a scaling to account for different epsilon
    
    std::cout << "Current cross section factor: " << this->GetCrossSecFactor() << std::endl;
}

void G4LDMEleScattering::SetAlphaD(G4double fac) {
// Set the factor to artificially increase the cross section

    alphaD = fac;
    std::cout << "The alphaD value for A' resonant production is set to " << alphaD << std::endl;
    
    double alphaD_0 = 0.5;
    double bias = alphaD/alphaD_0;
    std::cout << "Since the default value of alphaD is " << alphaD_0 << " , the cross section has to be scaled by a factor " << bias << std::endl;
    this->SetCrossSecFactor(CrossSecFactor*bias);
    
    std::cout << "Current cross section factor: " << this->GetCrossSecFactor() << std::endl;
}
