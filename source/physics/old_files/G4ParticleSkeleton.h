/*
 * G4Scalar.h
 *
 *  Created on: May 6, 2018
 *      Author: celentan
 */

#ifndef PHYSICS_G4SCALAR_H_
#define PHYSICS_G4SCALAR_H_

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                       G4Scalar                                 ###
// ######################################################################

class G4Scalar: public G4ParticleDefinition {
private:
	static G4Scalar* theInstance;
	G4Scalar() {
	}
	~G4Scalar() {
	}

	/*Since this is a singleton, the mass is used only at the first call*/
public:
	static G4Scalar* Definition(G4double mass=0., G4double g_mu=0);
	static G4Scalar* ScalarDefinition(G4double mass=0., G4double g_mu=0);
	static G4Scalar* Scalar(G4double mass=0., G4double g_mu=0);
};

#endif /* PHYSICS_G4SCALAR_H_ */
