/*
 * G4DarkMatter.h
 *
 *  Created on: March 15, 20220
 *      Author: celentan
 */

#ifndef PHYSICS_G4DarkMatter_H_
#define PHYSICS_G4DarkMatter_H_

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                       G4Scalar                                 ###
// ######################################################################

class G4DarkMatter: public G4ParticleDefinition {
private:
	static G4DarkMatter* theInstance;
	G4DarkMatter() {
	}
	~G4DarkMatter() {
	}

	/*Since this is a singleton, the mass is used only at the first call*/
public:
	static G4DarkMatter* Definition(G4double mass=0.03,G4int twoJ=2);
	static G4DarkMatter* DarkMatterDefinition(G4double mass=0.03,G4int twoJ=2);
	static G4DarkMatter* DarkMatter(G4double mass=0.03,G4int twoJ=2);
};

#endif /* PHYSICS_G4DarkMatter_H_ */
