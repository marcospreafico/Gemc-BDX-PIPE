/*
 * G4DarkMatter.cc
 *
 *  Created on: March 15, 2020
 *      Author: celentan
 */

#include "G4DarkMatter.h"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"



// ######################################################################
// ###                         G4DarkMatter                               ###
// ######################################################################

G4DarkMatter* G4DarkMatter::theInstance = 0;
G4DarkMatter* G4DarkMatter::Definition(G4double mass,G4int twoJ)
{
  if (theInstance !=0) return theInstance;
  const G4String name = "DarkMatter";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);



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
                 name,	mass, 	0,	 0,
				 twoJ,			1,		0,
				 0,			 0,                0,
				 "fermion",            0,        0,        51,
				 true,     0,             NULL,
				 false,      "fermion"
              );


  }
  theInstance = reinterpret_cast<G4DarkMatter*>(anInstance);
  return theInstance;
}

G4DarkMatter*  G4DarkMatter::DarkMatterDefinition(G4double mass,G4int twoJ)
{
  return Definition(mass,twoJ);
}

G4DarkMatter*  G4DarkMatter::DarkMatter(G4double mass,G4int twoJ)
{
  return Definition(mass,twoJ);
}
