/*
 * G4DarkPhotonProduction.cc
 *
 *  Created on: May 6, 2018
 *      Author: celentan
 */

#include "G4DarkPhotonAnnihilationProduction.h"

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
#include "G4Positron.hh"
#include "G4EmProcessSubType.hh"
#include "G4NistManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4DarkPhoton.h"
#include "G4DarkMatter.h"
#include "G4SystemOfUnits.hh"


G4DarkPhotonAnnihilationProduction::G4DarkPhotonAnnihilationProduction(const G4String &processName, G4ProcessType type) :
		G4VDiscreteProcess(processName, type), LowestEnergyLimit(MdarkPhoton), //the lower limit of the minimum energy
		HighestEnergyLimit(1e21 * eV), // ok to 1e21eV=1e12GeV, then LPM suppression
		CrossSecFactor(1.) {
	SetProcessSubType(G4DarkPhotonAnnihilationProductionProcessID);
	MeanFreePath = DBL_MAX;

	MdarkPhoton = (G4DarkPhoton::DarkPhoton()->GetPDGMass());
	WdarkPhoton = (G4DarkPhoton::DarkPhoton()->GetPDGWidth());
	MdarkMatter = (G4DarkMatter::DarkMatter()->GetPDGMass());

	TwoJdarkPhoton = G4DarkPhoton::DarkPhoton()->GetPDGiSpin();
	TwoJdarkMatter = G4DarkMatter::DarkMatter()->GetPDGiSpin();
	PdarkPhoton = G4DarkPhoton::DarkPhoton()->GetPDGiParity();

	//two "random" values to set from the external
	eps = 1E-4;
	alphaD = 0.5;



}

G4DarkPhotonAnnihilationProduction::~G4DarkPhotonAnnihilationProduction() {
}

G4bool G4DarkPhotonAnnihilationProduction::IsApplicable(const G4ParticleDefinition &particle) {
	return (&particle == G4Positron::Positron());
}

G4double G4DarkPhotonAnnihilationProduction::GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition*) {

// returns the particle mean free path in GEANT4 internal units
// (MeanFreePath is a private member of the class)


	const G4DynamicParticle *aDynamicParticle = aTrack.GetDynamicParticle();
	G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();
	G4Material *aMaterial = aTrack.GetMaterial();

	if (KineticEnergy <= LowestEnergyLimit) MeanFreePath = DBL_MAX;
	else
		MeanFreePath = ComputeMeanFreePath(KineticEnergy, aMaterial);

	return MeanFreePath;
}

G4double G4DarkPhotonAnnihilationProduction::ComputeMeanFreePath(G4double KineticEnergy, G4Material *aMaterial) {

// computes and returns the photon mean free path in GEANT4 internal units

	const G4ElementVector *theElementVector = aMaterial->GetElementVector();
	const G4double *NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

	G4double nSIGMA = 0;

	for (size_t i = 0; i < aMaterial->GetNumberOfElements(); ++i) {
		nSIGMA += NbOfAtomsPerVolume[i] * ComputeCrossSectionPerAtom(KineticEnergy, (*theElementVector)[i]);
	}
	return nSIGMA > DBL_MIN ? 1. / nSIGMA : DBL_MAX;
}

G4double G4DarkPhotonAnnihilationProduction::GetCrossSectionPerAtom(const G4DynamicParticle *aDynamicParticle, G4Element *anElement) {
	// gives the total cross section per atom in GEANT4 internal units
	G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();
	G4double crossSection = ComputeCrossSectionPerAtom(KineticEnergy, anElement);
	return crossSection;
}

//Ekin is in G4 internal units here!
G4double G4DarkPhotonAnnihilationProduction::ComputeCrossSectionPerAtom(G4double Ekin, G4Element *anElement) {

	// Calculates the microscopic cross section in GEANT4 internal units.
	// Note that this function also reads data "on the fly", storing data in the map

	if (Ekin <= eCut) return 0; // below threshold return 0
	G4double CrossSection = 0.0;
	G4double ZNucl = anElement->GetZ();

	//A.C. when doing the calculation, all the values are in G4 internal units
	G4double s = 2. * electron_mass_c2 * Ekin;
	if (sqrt(s) < 2. * MdarkMatter) return 0.;   // A.C. e+e- -> A' -> chi chi can happen also for an A' and chi with large mass,
	// i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed
	G4double q = sqrt(s) / 2. * sqrt(1 - 4 * MdarkMatter * MdarkMatter / (s));

	CrossSection = 4 * M_PI * fine_structure_const * eps * eps * alphaD;
	CrossSection = CrossSection * q / sqrt(s);
	CrossSection = CrossSection / ((s - MdarkPhoton * MdarkPhoton) * (s - MdarkPhoton * MdarkPhoton) + MdarkPhoton * MdarkPhoton * WdarkPhoton * WdarkPhoton);

	G4double K = 1;
	if (PdarkPhoton == -1) { //vector or pseudo-scalar
		if (TwoJdarkPhoton == 2) { //vector
			K = (s - 4. / 3. * q * q);

		} else if (TwoJdarkPhoton == 0) { //pseudo-scalar
			K = s / 2;
		}
	} else if (PdarkPhoton == 1) { //axial vector or scalar
		if (TwoJdarkPhoton == 2) { //axisl vector
			K = (8. / 3. * q * q);

		} else if (TwoJdarkPhoton == 0) { //scalar
			K = 2 * q * q;
		}
	}

	CrossSection = CrossSection * K;

	//hc1->Fill(Ekin/GeV);
	//hc2->Fill(Ekin/GeV,CrossSection/(1./GeV)/(1./GeV));



	//A.C. correct here for atomic effects
	CrossSection = CrossSection * ZNucl;

	// increase the CrossSection by CrossSecFactor (by default 1)
	CrossSection *= CrossSecFactor;

	// units: the dimensions of CrossSection up to this point are inverse energy squared. I need to multiply by (hbar*c)^2 to have it as a length^2
	CrossSection = CrossSection * hbarc_squared;
	return CrossSection;
}

G4VParticleChange* G4DarkPhotonAnnihilationProduction::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {

	aParticleChange.Initialize(aTrack);
	G4Material *aMaterial = aTrack.GetMaterial();

// current muon kinetic energy and direction, return if energy too low
	const G4DynamicParticle *aDynamicParticle = aTrack.GetDynamicParticle();
	G4double kineticE = aDynamicParticle->GetKineticEnergy();
	if (kineticE <= LowestEnergyLimit) {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

// select randomly one element constituting the material
	const G4Element *anElement = SelectRandomAtom(aDynamicParticle, aMaterial);
//try again
	if (anElement == NULL) {
		anElement = SelectRandomAtom(aDynamicParticle, aMaterial);
	}
//if still bad, do nothing
	if (anElement == NULL) {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	/*Define the A' momentum as the dark photon momentum*/
	G4ThreeVector ImpingingMomentum = aDynamicParticle->GetMomentum();
	G4ThreeVector DarkPhotonDirection=ImpingingMomentum.unit();


	aParticleChange.SetNumberOfSecondaries(1);
// create G4DynamicParticle object for the particle1
	G4DynamicParticle *aParticle1 = new G4DynamicParticle(G4DarkPhoton::DarkPhoton(), DarkPhotonDirection, kineticE);
	aParticleChange.AddSecondary(aParticle1);

	aParticleChange.ProposeTrackStatus(fStopAndKill);

	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4Element* G4DarkPhotonAnnihilationProduction::SelectRandomAtom(const G4DynamicParticle *aDynamicParticle, G4Material *aMaterial) {
// select randomly 1 element within the material, invoked by PostStepDoIt

	const G4int NumberOfElements = aMaterial->GetNumberOfElements();
	const G4ElementVector *theElementVector = aMaterial->GetElementVector();
	if (NumberOfElements == 1) return (*theElementVector)[0];

	const G4double *NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

	G4double PartialSumSigma = 0.;
	G4double rval = G4UniformRand() / MeanFreePath;

	for (G4int i = 0; i < NumberOfElements; ++i) {
		PartialSumSigma += NbOfAtomsPerVolume[i] * GetCrossSectionPerAtom(aDynamicParticle, (*theElementVector)[i]);
			if (rval <= PartialSumSigma) return ((*theElementVector)[i]);
	}
	std::cout << " WARNING !!! - The Material '" << aMaterial->GetName() << "' has no elements, NULL pointer returned. rval: " << rval << std::endl;
	return NULL;
}

void G4DarkPhotonAnnihilationProduction::SetCrossSecFactor(G4double fac) {
// Set the factor to artificially increase the cross section

	CrossSecFactor = fac;
	std::cout << "The cross section for A' resonant production is artificially " << "increased by the CrossSecFactor=" << CrossSecFactor << std::endl;
}

void G4DarkPhotonAnnihilationProduction::SetEps(G4double fac) {
// Set the factor to artificially increase the cross section

	eps = fac;
	std::cout << "The eps value for A' resonant production is set to " << eps << std::endl;
}

void G4DarkPhotonAnnihilationProduction::SetAlphaD(G4double fac) {
// Set the factor to artificially increase the cross section

	alphaD = fac;
	std::cout << "The alphaD value for A' resonant production is set to " << alphaD << std::endl;
}

void G4DarkPhotonAnnihilationProduction::SetEnergyCut(G4double fac) {
// Set the factor to artificially increase the cross section
	eCut = fac;
	std::cout << "The energy value for A' resonant production is set to " << eCut/GeV <<" GeV "<< std::endl;
}
