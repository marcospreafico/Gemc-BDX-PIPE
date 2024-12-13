/*
 * G4DarkPhoton.h
 *
 *  Created on: March 15, 20220
 *      Author: celentan
 */

#ifndef PHYSICS_G4DARKPHOTON_H_
#define PHYSICS_G4DARKPHOTON_H_

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                       G4Scalar                                 ###
// ######################################################################

class G4DarkPhoton: public G4ParticleDefinition {
private:
	static G4DarkPhoton* theInstance;
	G4DarkPhoton() {
	}
	~G4DarkPhoton() {
	}

	/*Since this is a singleton, the mass is used only at the first call*/
public:
	static G4DarkPhoton* Definition(G4double mass=0.1,G4int twoJ=2,G4int P=-1,G4double eps=1E-4,G4double alphaD=0.1);
	static G4DarkPhoton* DarkPhotonDefinition(G4double mass=0.,G4int twoJ=2,G4int P=-1,G4double eps1=1E-4,G4double alphaD=0.1);
	static G4DarkPhoton* DarkPhoton(G4double mass=0.1,G4int twoJ=2,G4int P=-1,G4double eps=1E-4,G4double alphaD=0.1);
};

#endif /* PHYSICS_G4DARKPHOTON_H_ */
