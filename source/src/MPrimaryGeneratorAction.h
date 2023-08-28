#ifndef MPrimaryGeneratorAction_h
#define MPrimaryGeneratorAction_h 1

// Geant4
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"

// For reading StdHep files
#include "lStdHep.hh"
using namespace UTIL;

// gemc
#include "options.h"

// C++
#include <fstream>

class TFile;
class TH3D;
class TH2D;
class TTree;

class MEventAction;

class userInforForParticle {
public:
	vector<double> infos;
};

// Class definition
class MPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
	MPrimaryGeneratorAction(goptions*,MEventAction *evtAction);
	~MPrimaryGeneratorAction();

public:
	void GeneratePrimaries(G4Event* anEvent);
	goptions *gemcOpt;
	double GEN_VERBOSITY;
	double getBeamPol() {
		return beamPol;
	}

	// user defined info for header
	vector<double> headerUserDefined;   ///< user defined infos in the  header

	// user defined info for each particle
	vector<userInforForParticle> userInfo;

	// these two should belong to a general generator class
	// start time can be modeled
	double getTimeWindow() {
		return TWINDOW;
	}
	double getStartTime() {
		return TWINDOW / 2;
	}

	bool isFileOpen() {
		return !gif.eof();
	}

	void setNsup(int m_nsup) {
		Nsup = m_nsup;
	}
	int getNsup() {
		return Nsup;
	}

private:
	string input_gen;                 ///< Input Option: internal or external
	string background_gen;            ///< Input Option: background from file in LUND format
	string cosmics;                   ///< cosmic ray option
	string cArea;
	string JLabFlux;				  ///< JLabFlux - muons sampling option
	string JLabFlux2;				  ///< JLabFlux - muons sampling option
	string hd_msg;                    ///< Head Message Log
	int ntoskip;                      ///< Number of events to skip
	int eventIndex;                   ///< Set to 1

	G4ParticleTable* particleTable;   ///< Geant4 Particle Table

	// Primary Beam
	G4ParticleDefinition *Particle;   ///< Particle type
	double mom, dmom;               ///< beam momentum, delta momentum
	double theta, dtheta;             ///< theta, delta theta
	double phi, dphi;               ///< phi, delta phi
	double vx, vy, vz;                ///< Beam Vertex coordinates
	double dvr, dvz;                  ///< Deltas Beam Vertex: Radius and z-vertex
	double dvx, dvy;                 ///< Deltas Beam Vertex: individual coordinates
	int gaussOrFlatV;                 ///< 0 (default): flat distribution 1: gaussian distribution
	int drdzOrdxdydz;                 ///< spread is in (dr, dz): 0. (dx, dy dz): 1
	double polDeg, polTheta, polPhi;  ///< Polarization degree and  direction
	double ctheta;                    ///< customized theta direction for the z axis
	double cphi;                      ///< customize phi direction for the z axis
	double primaryFlat;               ///< if this is set to 1 spread flat in theta, not cos(theta)

	double cosmicA, cosmicB, cosmicC; ///< cosmic ray model parameters
	double cminp, cmaxp, cMom;        ///< minimum and maximum cosmic ray momentum
	G4ThreeVector cosmicTarget;       ///< Location of area of interest for cosmic rays
	double cosmicRadius;              ///< radius of area of interest for cosmic rays
	double cosmicRadiusBackProj;      ///< radius of back-projection
	G4ThreeVector HallDim;              ///< Hall diemension threev ector
	double HallRadius;              ///< Hall radius
	string cosmicGeo;                 ///< type of surface for cosmic ray generation (sphere || cylinder)
	string cosmicParticle;            ///< type of cosmic ray particle (muon || neutron)
	int muonDecay;                  ///< type of muon decay

	//JLabFlux
	int JLabFluxID;
	double JLabFluxDetHeight;
	double JLabFluxMinE;
	string JLabFluxROOTfname;
	TFile *JLabFluxROOTfile;
	TH3D *JLabFluxROOThisto;
	TTree *JLabFluxROOTtree;

	double cxMu, cyMu, czMu; //cosine directors;
	double kinMu; //kinetic energy in GeV
	double vxMu, vyMu, vzMu; //vertex in cm
	double wMu; //weight
	int idFluxMu;
	int idxMu;

	double zOrigin;
	double zPipe1, zPipe2;
	double JLabFluxDeltaZ;
	double JLabFluxDeltaX;
	double JLabFluxDeltaY;
	double JLabFluxNORM;
	double JLabFluxEFF;

	// Generators Input Files
	ifstream gif;                    ///< Generator Input File
	ifstream bgif;                   ///< Background Generator Input File
	string gformat;                ///< Generator Format. Supported: LUND.
	string gfilename;              ///< Input Filename for main events
	double beamPol;                ///< Beam Polarization as from the LUND format, it
	int beagleHeader;                 ///< Starts at 0, will become 1 once it is read

	lStdHep *stdhep_reader;         /// Handle to the object for reading StdHep files.

	// Luminosity Beam
	G4ParticleDefinition *L_Particle;  ///< Luminosity Particle type
	double L_mom, L_dmom;             ///< Luminosity beam momentum, delta momentum
	double L_theta, L_dtheta;          ///< Luminosity theta,  delta theta
	double L_phi, L_dphi;              ///< Luminosity phi, delta phi, randomized phi
	double L_vx, L_vy, L_vz;           ///< Luminosity Beam Vertex coordinates
	double L_dvr, L_dvz;              ///< Luminosity Deltas Beam Vertex: Radius and z-vertex
	int NP;                            ///< Number of Luminosity Particles per event
	double TWINDOW;                    ///< Time Window
	double TBUNCH;                     ///< Time Between Bunches
	double lumiFlat;                   ///< if this is set to 1, spread flat in theta, not cos(theta)

	// Luminosity Beam2
	G4ParticleDefinition *L2_Particle;    ///< Luminosity Particle type
	double L2_mom, L2_dmom;               ///< Luminosity beam momentum, delta momentum
	double L2_theta, L2_dtheta;           ///< Luminosity theta,  delta theta
	double L2_phi, L2_dphi;               ///< Luminosity phi, delta phi, randomized phi
	double L2_vx, L2_vy, L2_vz;           ///< Luminosity Beam Vertex coordinates
	double L2_dvr, L2_dvz;                ///< Luminosity Deltas Beam Vertex: Radius and z-vertex
	int NP2;                              ///< Number of Luminosity Particles per event
	double TBUNCH2;                       ///< Time Between Bunches
	double lumi2Flat;                     ///< if this is set to 1, spread flat in theta, not cos(theta)


	G4ParticleGun* particleGun;
	void setBeam();

	double cosmicMuBeam(double, double);
	double cosmicNeutBeam(double, double);

	int Nsup;
	int Ntot;

	void setParticleFromPars(int, int, int, int, double, double, double, double, double, double, G4Event* anEvent);

	MEventAction* evtAction;

};

#endif

