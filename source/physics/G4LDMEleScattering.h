/*
 * G4LDMEleScattering.h
 *
 *  Created on: March 18th, 2024
 *      Author: mspreafico
 */

#ifndef PHYSICS_G4LDMEleScattering_H_
#define PHYSICS_G4LDMEleScattering_H_

#define G4LDMEleScatteringProcessID 800 //This is the process id to be used in SetProcessSubType. Is it used?

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Element.hh"
#include "G4Electron.hh" // !! Need to include the recoiling particle
#include "G4Step.hh"

class TGraph;
class TFile;
class TH2D;
class TObject;
class TKey;
class G4LDMEleScattering: public G4VDiscreteProcess {
public:
	G4LDMEleScattering(const G4String& processName = "LDMEleScattering", G4ProcessType type = fElectromagnetic);
	virtual ~G4LDMEleScattering();

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
	//G4double ComputeCrossSectionPerAtom(G4double KineticEnergy, G4String ElementSymbol);
    G4double ComputeCrossSectionPerAtom(G4double KineticEnergy, G4Element* anElement);
    
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

private:

	G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicParticle, G4Material* aMaterial);
	TFile* GetDataFile();

private:

	G4double MDM;                   // Dark Matter mass 
	G4double LowestEnergyLimit;     // low  energy limit of the tables
	G4double HighestEnergyLimit;    // high energy limit of the tables

	G4double MeanFreePath;           // actual MeanFreePath (current medium)
	G4double CrossSecFactor;         // factor to artificially increase

	//The path to the dataFolder ($BDX_SCALAR_DATA/mu)
	G4String dataPathFolder;

    G4double eps = 1e-3; 
    G4double alphaD; 
    
	TFile *fData;
	TGraph *gXsecData;
	TH2D *hXsecDistr;

};

#endif /* PHYSICS_G4LDMEleScattering_H_ */
