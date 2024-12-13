/*
 * G4ScalarProduction.cc
 *
 *  Created on: May 6, 2018
 *      Author: celentan
 */

#include "G4MuonScalarProduction.h"

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

#include "G4Scalar.h"
#include "G4SystemOfUnits.hh"

//root headers
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
#include "TObject.h"
#include "TList.h"
#include "TKey.h"

//Keeping global variables shared possibly by different istances of the same file means this is not thread-safe.
//Will use gemc in single-thread mode.

static bool messageOut = true;

//This holds the total cross section, as a function of the impinging muon KINETIC energy, for the given element.
//I can store data as any "table". I use root TGraph since it has built-in interpolation features
static std::map<G4String, TGraph *> xsecValData;

G4MuonScalarProduction::G4MuonScalarProduction(const G4String& processName, G4ProcessType type) :
		G4VDiscreteProcess(processName, type), Mscalar(G4Scalar::Scalar()->GetPDGMass()), LowestEnergyLimit(Mscalar), //the lower limit of the minimum energy
		HighestEnergyLimit(1e21 * eV), // ok to 1e21eV=1e12GeV, then LPM suppression
		CrossSecFactor(1.), fData(0), gXsecData(0), hXsecDistr(0) {
	SetProcessSubType(G4MuonScalarProductionProcessID);
	MeanFreePath = DBL_MAX;

	/*The data is loaded "on the fly", in the GetMeanFreePath - for the given ELEMENT.*/
	char* pPath;
	pPath = getenv("BDX_SCALAR_DATA");
	if (pPath == NULL) {
		std::cout << std::endl << std::endl;
		std::cout << "G4MuonScalarProduction creator, the environmental variable BDX_SCALAR_DATA is not set." << std::endl;
		std::cout << "This should be set to the folder with the data for the scalar production, organized as a set of folders as follows:" << std::endl;
		std::cout << "First level: particle name (muon/electron). In this case, muon is relevant" << std::endl;
		std::cout << "Second level: the ROOT files with the data, as data_ElementSymbol_ScalarMassInGeV.root" << std::endl;
		std::cout << "The scalarMassInGeV should have a precision of 3 (i.e. 0.010 is ok, 0.01 is not for 10 MeV)" << std::endl;
		exit(1);
	} else {
		dataPathFolder = G4String(pPath);
		std::cout << std::endl << std::endl;
		if (messageOut) {
			std::cout << "G4MuonScalarProduction creator, the environmental variable BDX_SCALAR_DATA is set to: " << dataPathFolder << std::endl;
			messageOut = false;
		}
		dataPathFolder += "/muon";
		struct stat sb;
		if (!(stat(dataPathFolder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) {
			std::cout << std::endl << std::endl;
			std::cout << "G4MuonScalarProduction creator, muon folder is missing within BDX_SCALAR_DATA" << std::endl;
			exit(1);
		}

	}
}

G4MuonScalarProduction::~G4MuonScalarProduction() {
	// TODO Auto-generated destructor stub
}

G4bool G4MuonScalarProduction::IsApplicable(const G4ParticleDefinition& particle) {
	return ((&particle == G4MuonPlus::MuonPlus()) || ((&particle == G4MuonMinus::MuonMinus())));
}

G4double G4MuonScalarProduction::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*) {

// returns the particle mean free path in GEANT4 internal units
// (MeanFreePath is a private member of the class)

	aTrack.GetParticleDefinition()->GetPDGMass();

	const G4DynamicParticle* aDynamicParticle = aTrack.GetDynamicParticle();
	G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();
	G4double Mass = aTrack.GetParticleDefinition()->GetPDGMass();
	G4Material* aMaterial = aTrack.GetMaterial();

	if (KineticEnergy <= LowestEnergyLimit) MeanFreePath = DBL_MAX;
	else
		MeanFreePath = ComputeMeanFreePath(KineticEnergy, aMaterial);

	return MeanFreePath;
}

G4double G4MuonScalarProduction::ComputeMeanFreePath(G4double KineticEnergy, G4Material* aMaterial) {

// computes and returns the photon mean free path in GEANT4 internal units

	const G4ElementVector* theElementVector = aMaterial->GetElementVector();
	const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

	G4double nSIGMA = 0;

	for (size_t i = 0; i < aMaterial->GetNumberOfElements(); ++i) {
		G4double AtomicZ = (*theElementVector)[i]->GetZ();
		G4double AtomicA = (*theElementVector)[i]->GetA() / (g / mole);
		G4String AtomicSymbol = (*theElementVector)[i]->GetSymbol();
		nSIGMA += NbOfAtomsPerVolume[i] * ComputeCrossSectionPerAtom(KineticEnergy, AtomicSymbol);
	}
	return nSIGMA > DBL_MIN ? 1. / nSIGMA : DBL_MAX;
}

G4double G4MuonScalarProduction::GetCrossSectionPerAtom(const G4DynamicParticle* aDynamicParticle, G4Element* anElement)

// gives the total cross section per atom in GEANT4 internal units
		{
	G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();
	G4String AtomicSymbol = anElement->GetSymbol();
	G4double crossSection = ComputeCrossSectionPerAtom(KineticEnergy, AtomicSymbol);
	return crossSection;
}

G4double G4MuonScalarProduction::ComputeCrossSectionPerAtom(G4double Ekin, G4String ElementSymbol) {

	// Calculates the microscopic cross section in GEANT4 internal units.
	// Note that this function also reads data "on the fly", storing data in the map

	if (Ekin <= LowestEnergyLimit) return 0; // below threshold return 0

	std::map<G4String, TGraph *>::iterator it;
	it = xsecValData.find(ElementSymbol);

	if (it == xsecValData.end()) {
		//This means the material was not found already. Search for it.

		//std::cout << "G4MuonScalarProduction::ComputeCrossSectionPerAtom, element: " << ElementSymbol << " not in the cache. Search it" << std::endl;
	  fData = this->GetDataFile(ElementSymbol);
		
	  if (fData->IsZombie()) {
		  //std::cout << "File does not exists! Return 0 xsection" << std::endl;
			delete fData;
			return 0.0;
		}
		
		gXsecData = NULL;
		gXsecData = (TGraph*) fData->Get("gData");
		if (gXsecData == NULL) {
			//std::cout << "Data does not exists! Return 0 xsection" << std::endl;
			delete fData;
			return 0.0;
		}
		/*I clone the original data graph, then delete the file*/
		xsecValData[ElementSymbol] = (TGraph*) gXsecData->Clone("gData");
		delete fData;
	} else {
		gXsecData = it->second;
	}

	G4double CrossSection = 0.0;

	CrossSection = gXsecData->Eval(Ekin/GeV); //Need to evaluate in GeV!!!
	CrossSection = CrossSection * picobarn; //data is in pbarn
	//std :: cout << CrossSecFactor << std::endl;
	CrossSection *= CrossSecFactor; // increase the CrossSection by  (by default 1)
	return CrossSection;
}

G4VParticleChange* G4MuonScalarProduction::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
//
// generation of scalar production
//
  //std::cout << "Scalar production" << std::endl;

	G4double deltaMin = DBL_MAX;
	G4double delta;
	G4int idx;

	G4double randomE, randomZ, randomPhi;

	aParticleChange.Initialize(aTrack);
	G4Material* aMaterial = aTrack.GetMaterial();

	// current muon kinetic energy and direction, return if energy too low
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


	//Get the atomic symbol
	G4String ElementSymbol = anElement->GetSymbol();

	std::map<G4String, TGraph *>::iterator it;
	it = xsecValData.find(ElementSymbol);

	if (it == xsecValData.end()) {
		std::cerr << "Error in G4MuonScalarProduction::PostStepDoIt(), element " << ElementSymbol << " not found. Doing nothing " << std::endl;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	//Get the graph with the data
	gXsecData = it->second;
	//Compute the index for this kinetic energy Ekin -> the index of the closest point. I assume nothing about the points (not equally spaced, nor ordered)
	for (int ip = 0; ip < gXsecData->GetN(); ip++) {
		delta = fabs(kineticE / GeV - gXsecData->GetX()[ip]);
		if (delta < deltaMin) {
			deltaMin = delta;
			idx = ip;
		}
	}

	/*Get the fData and from this get the histogram we need*/
	fData = this->GetDataFile(ElementSymbol);
	if (fData->IsZombie()) {
		//std::cout << "File does not exists! Return 0 xsection" << std::endl;
		delete fData;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}
	hXsecDistr = (TH2D*) fData->Get(Form("h2_%i", idx));

	//Following error may be due to binning+threshold effects
	if (hXsecDistr == NULL) {
		//std::cout << "Error! Element " << ElementSymbol << " not finding histo. Doing nothing. Ekin: " << kineticE << " idx: " << idx << std::endl;
		delete fData;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	/*If we are here, it means we have the histogram. Get random number from it: kinetic energy in GeV (X) and cosine of scattering angle (Y)*/
	hXsecDistr->GetRandom2(randomE, randomZ);
	int repeat_extraction = 0; 
	while((randomE*GeV)+Mscalar > kineticE && repeat_extraction<10){
	  hXsecDistr->GetRandom2(randomE, randomZ);
	  std::cout << "Exceeding muon energy, repeat extraction" << std::endl;
	  repeat_extraction ++; 
	}
	randomE = randomE * GeV; //always need to specify the units in G4.
	if(repeat_extraction == 9) {
	  std::cout << "End extraction, could not find correct combination" << std::endl; 
	}

	//Just in case..
	if (randomZ > 1) randomZ = 1;
	if (randomZ < -1) randomZ = -1;
	randomPhi = G4UniformRand()*twopi;
	G4ThreeVector ScalarDirection(sqrt(1 - randomZ * randomZ) * sin(randomPhi), sqrt(1 - randomZ * randomZ) * cos(randomPhi), randomZ);
	// rotate to actual impinging muon direction
	ScalarDirection.rotateUz(MomentumDirection);

	/*Define the scattered muon momentum as the difference: impinging - scalar*/
	G4ThreeVector ImpingingMomentum = aDynamicParticle->GetMomentum();
	G4ThreeVector ScalarMomentum = ScalarDirection * sqrt(randomE * randomE + 2 * randomE * Mscalar);
	G4ThreeVector MuonMomentum = ImpingingMomentum - ScalarMomentum;

	std::cout << randomE << " " << GeV << " " << randomE/GeV << " " << Mscalar << std::endl;

	std::cout << "random extraction done, impinging energy: " << kineticE / GeV << " kinetic energy " << randomE / GeV << " dir: " << randomZ << " mass " << Mscalar << std::endl;
	//std::cout << ImpingingMomentum.x() << " " << ImpingingMomentum.y() << " " << ImpingingMomentum.z() << std::endl;
	std::cout << ScalarMomentum.x() << " " << ScalarMomentum.y() << " " << ScalarMomentum.z() << std::endl;
	//std::cout << MuonMomentum.x() << " " << MuonMomentum.y() << " " << MuonMomentum.z() << std::endl;

	aParticleChange.SetNumberOfSecondaries(1);
	// create G4DynamicParticle object for the particle1
	G4DynamicParticle* aParticle1 = new G4DynamicParticle(G4Scalar::Scalar(), ScalarDirection, randomE);
	aParticleChange.AddSecondary(aParticle1);

	
	const G4ThreeVector* VertexPosition =  aParticleChange.GetPosition();
	std::cout << VertexPosition->x() << " " << VertexPosition->y() << " " << VertexPosition->z() << std::endl;
	const G4double weight = aDynamicParticle->GetPrimaryParticle()->GetTotalEnergy(); 
	std::cout << weight << " <-- primary E" << std::endl;

	// Change the incident particle
	aParticleChange.ProposeMomentumDirection(MuonMomentum.unit());
	G4double finalKE = (kineticE -((randomE)+Mscalar)); 
	if(finalKE<0){
	  std::cout << "WARNING: negative final kinetic energy - " << finalKE << std:: endl;
	  finalKE = 0; 
	}
	aParticleChange.ProposeEnergy(kineticE -((randomE)+Mscalar));
	std::cout << kineticE << std::endl;

	//aParticleChange.ProposeTrackStatus(fStopAndKill);

	//  Reset NbOfInteractionLengthLeft and return aParticleChange
	delete fData;


	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4Element* G4MuonScalarProduction::SelectRandomAtom(const G4DynamicParticle* aDynamicParticle, G4Material* aMaterial) {
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

void G4MuonScalarProduction::SetCrossSecFactor(G4double fac) {
// Set the factor to artificially increase the cross section

	CrossSecFactor = fac;
	//std::cout << "The cross section for ScalarProductin is artificially " << "increased by the CrossSecFactor=" << CrossSecFactor << std::endl;
}

/*This retuns a new TFile for the current element symbol.
 * It is up to the user to delete it!
 */
TFile* G4MuonScalarProduction::GetDataFile(G4String ElementSymbol) {
	G4double MscalarTMP = Mscalar / GeV;
	G4String massString;
	G4String filePath;

	std::stringstream stream;
	stream << std::fixed << std::setprecision(3) << MscalarTMP;

	massString = G4String(stream.str().c_str());
	filePath = dataPathFolder + "/data_" + ElementSymbol + "_" + massString + ".root";

	return new TFile(filePath.c_str());
}
