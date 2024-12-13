/*
 * G4DarkPhotonAnnihilationProduction.cc
 *
 *  Created on: Mar 15, 2021
 *      Author: celentan
 */

#ifndef PHYSICS_G4DARKPHOTONANNIHILATIONPRODUCTION_CC_
#define PHYSICS_G4DARKPHOTONANNIHILATIONPRODUCTION_CC_

#define G4DarkPhotonAnnihilationProductionProcessID 101 //This is the process id to be used in SetProcessSubType. Is it used?

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Element.hh"
#include "G4Positron.hh"
#include "G4Step.hh"


#include "TH1D.h"
#include "TFile.h"


class G4DarkPhotonAnnihilationProduction: public G4VDiscreteProcess {
public:
	G4DarkPhotonAnnihilationProduction(const G4String &processName = "DarkPhotonAnnihilationProduction", G4ProcessType type = fElectromagnetic);
	virtual ~G4DarkPhotonAnnihilationProduction();

	G4bool IsApplicable(const G4ParticleDefinition&);
	// true for e+ for the moment

	G4double GetMeanFreePath(const G4Track &aTrack, G4double previousStepSize, G4ForceCondition *condition);
	// It returns the MeanFreePath of the process for the current track :
	// (energy, material)
	// The previousStepSize and G4ForceCondition* are not used.
	// This function overloads a virtual function of the base class.
	// It is invoked by the ProcessManager of the Particle.

	G4double GetCrossSectionPerAtom(const G4DynamicParticle *aDynamicGamma, G4Element *anElement);
	// It returns the total CrossSectionPerAtom of the process,
	// for the current DynamicGamma (energy), in anElement.

	G4VParticleChange* PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);
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

	void SetEps(G4double fac);
	// Set the factor to artificially increase the crossSection (default 1)

	G4double GetEps() {
		return eps;
	}
	// Get the factor to artificially increase the cross section

	void SetAlphaD(G4double fac);
	// Set the factor to artificially increase the crossSection (default 1)

	G4double GetAlphaD() {
		return alphaD;
	}

	G4double GetEnergyCut(){
	  return eCut;
	}
	
	void SetEnergyCut(G4double fac);

	// Get the factor to artificially increase the cross section

	G4double ComputeMeanFreePath(G4double KineticEnergy, G4Material *aMaterial);
	G4double ComputeCrossSectionPerAtom(G4double KineticEnergy, G4Element *anElement);

private:

	G4Element* SelectRandomAtom(const G4DynamicParticle *aDynamicParticle, G4Material *aMaterial);

private:

	G4double MdarkPhoton;
	G4double WdarkPhoton;
	G4double MdarkMatter;
	G4double eps;
	G4double alphaD;
	G4int TwoJdarkPhoton;
	G4int PdarkPhoton;
	G4int TwoJdarkMatter;

	G4double LowestEnergyLimit;     // low  energy limit of the tables
	G4double HighestEnergyLimit;    // high energy limit of the tables

	G4double MeanFreePath;           // actual MeanFreePath (current medium)
	G4double CrossSecFactor;         // factor to artificially increase

	G4double eCut;


};

#endif /* PHYSICS_G4DARKPHOTONANNIHILATIONPRODUCTION_CC_ */
