/*
 * G4ScalarProduction.h
 *
 *  Created on: May 6, 2018
 *      Author: celentan
 */

#ifndef PHYSICS_G4MUONSCALARPRODUCTION_H_
#define PHYSICS_G4MUONSCALARPRODUCTION_H_

#define G4MuonScalarProductionProcessID 100 //This is the process id to be used in SetProcessSubType. Is it used?

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Element.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Step.hh"

class TGraph;
class TFile;
class TH2D;
class TObject;
class TKey;
class G4MuonScalarProduction: public G4VDiscreteProcess {
public:
	G4MuonScalarProduction(const G4String& processName = "ScalarProduction", G4ProcessType type = fElectromagnetic);
	virtual ~G4MuonScalarProduction();

	G4bool IsApplicable(const G4ParticleDefinition&);
	// true for muon for the moment

	G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition);
	// It returns the MeanFreePath of the process for the current track :
	// (energy, material)
	// The previousStepSize and G4ForceCondition* are not used.
	// This function overloads a virtual function of the base class.
	// It is invoked by the ProcessManager of the Particle.

	G4double GetCrossSectionPerAtom(const G4DynamicParticle* aDynamicGamma, G4Element* anElement);
	// It returns the total CrossSectionPerAtom of the process,
	// for the current DynamicGamma (energy), in anElement.

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
	// It computes the final state of the process (at end of step),
	// returned as a ParticleChange object.
	// This function overloads a virtual function of the base class.
	// It is invoked by the ProcessManager of the Particle.

	void SetCrossSecFactor(G4double fac);
	// Set the factor to artificially increase the crossSection (default 1)

	G4double GetCrossSecFactor() {
		return CrossSecFactor;
	}
	// Get the factor to artificially increase the cross section

	G4double ComputeMeanFreePath(G4double KineticEnergy, G4Material* aMaterial);
	G4double ComputeCrossSectionPerAtom(G4double KineticEnergy, G4String ElementSymbol);

private:

	G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicParticle, G4Material* aMaterial);
	TFile* GetDataFile(G4String ElementSymbol);

private:

	G4double Mscalar;
	G4double LowestEnergyLimit;     // low  energy limit of the tables
	G4double HighestEnergyLimit;    // high energy limit of the tables

	G4double MeanFreePath;           // actual MeanFreePath (current medium)
	G4double CrossSecFactor;         // factor to artificially increase

	//The path to the dataFolder ($BDX_SCALAR_DATA/mu)
	G4String dataPathFolder;

	TFile *fData;
	TGraph *gXsecData;
	TH2D *hXsecDistr;

};

#endif /* PHYSICS_G4MUONSCALARPRODUCTION_H_ */
