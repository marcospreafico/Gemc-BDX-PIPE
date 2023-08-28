/*
 * G4DarkPhoton.cc
 *
 *  Created on: March 15, 2020
 *      Author: celentan
 */

#include "G4DarkPhoton.h"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                         G4DarkPhoton                               ###
// ######################################################################

G4DarkPhoton* G4DarkPhoton::theInstance = 0;
G4DarkPhoton* G4DarkPhoton::Definition(G4double mass,G4int twoJ,G4int P,G4double eps,G4double alphaD)
{
	if (theInstance != 0) return theInstance;
	// search in particle table
	G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();

	//search the DM to get the mass
	const G4String DMname = "DarkMatter";
	G4ParticleDefinition *anInstanceDM = pTable->FindParticle(DMname);
	if (anInstanceDM == 0) {
		G4cerr << " ERROR in G4DarkPhoton::Definition, DarkMatter not found. You should construct it first" << G4endl;
		exit(1);
	}
	G4double DMmass = anInstanceDM->GetPDGMass();
	G4int DMJ = anInstanceDM->GetPDGiSpin(); //pay attention, GetPDGiSpin() returns the spin in units of 1/2;

	G4double ratio = DMmass / mass;
	const G4String name = "DarkPhoton";
	G4ParticleDefinition *anInstance = pTable->FindParticle(name);

	G4double width, lifetime;
	if (DMJ == 1) {
		if (P == -1) { //vector and pseudo-scalar
			if (twoJ == 2) { //vector
				width = alphaD * mass / 3 * sqrt(1 - 4 * ratio * ratio) * (1 + 2 * ratio * ratio);
			}
			else if (twoJ == 0){ //pseudo-scalar
				width = alphaD * mass / 2 * sqrt(1 - 4 * ratio * ratio);
			}
			else{
				G4cerr << " DP 2*J is: " << twoJ << " parity negative: not yet supported " << G4endl;
				exit(1);
			}
		}
		else if (P == 1) { //axial vector and scalar
			if (twoJ == 2) { //axial vector
				width = alphaD * mass / 3 * pow(1 - 4 * ratio * ratio,3./2);
			} else if (twoJ == 0) { //scalar
				width = alphaD * mass / 2 * pow(1 - 4 * ratio * ratio,3./2);
			} else {
				G4cerr << " DP 2*J is: " << twoJ << " parity positive: not yet supported " << G4endl;
				exit(1);
			}
		}else{
			G4cerr << "DP Parity: " << P << " not supported" << G4endl;
			exit(1);
		}
	}else{
		G4cerr << " DarkMatter 2*J is: " << DMJ << " not supported" << G4endl;
		exit(1);
	}

	lifetime=(hbar_Planck/width); //in time units


  if (anInstance ==0)
  {
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  anInstance = new G4ParticleDefinition(
                 name,	mass, 	width,	 0,
				 twoJ,     P,		0,
				 0,			 0,                0,
				 "boson",            0,        0,        55,
				 false,     lifetime,             NULL,
				 false,      "vector"
              );

  //create Decay Table
  G4DecayTable* table = new G4DecayTable();
  // create decay channels
  // A' -> DM DM
  G4VDecayChannel* mode  = new G4PhaseSpaceDecayChannel("DarkPhoton",1.000,2,"DarkMatter","DarkMatter");
  table->Insert(mode);
  anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<G4DarkPhoton*>(anInstance);
  return theInstance;
}

G4DarkPhoton*  G4DarkPhoton::DarkPhotonDefinition(G4double mass,G4int twoJ,G4int P,G4double eps,G4double alphaD)
{
  return Definition(mass,twoJ,P,eps,alphaD);
}

G4DarkPhoton*  G4DarkPhoton::DarkPhoton(G4double mass,G4int twoJ,G4int P,G4double eps,G4double alphaD)
{
  return Definition(mass,twoJ,P,eps,alphaD);
}
