// c++ headers
#include <iostream>

#include <cstdlib>

#include <vector>

using namespace std;

// gemc headers
#include "PhysicsList.h"
#include "PhysicsListMessenger.h"
#include "string_utilities.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// geant4 headers
#include "G4LossTableManager.hh"
#include "G4PhysListFactory.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4DecayTable.hh"
#include "G4ProcessTable.hh"

// geant4 physics headers
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4GammaConversionToMuons.hh"

#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT_TRV.hh"
#include "G4HadronPhysicsFTF_BIC.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4HadronPhysicsShielding.hh"


//DarkPhoton
#include "G4DarkPhotonAnnihilationProduction.h"
#include "G4DarkMatter.h"
#include "G4DarkPhoton.h"

#include "G4HadronicProcessStore.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4ProcessVector.hh"
#include "G4HadronicProcessType.hh"


// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

PhysicsList::PhysicsList(goptions opts) :
		G4VModularPhysicsList() {


	gemcOpt = opts;
	verbosity = gemcOpt.optMap["PHYS_VERBOSITY"].arg;
	ingredientsList = gemcOpt.optMap["PHYSICS"].args;
	int fastmcMode = gemcOpt.optMap["FASTMCMODE"].arg;

	// default physics lists
	hadronicPhys = "none";
	EMPhys = "none";
	opticalPhys = "none";
	HPPhys = "none";

	// loading hadronic and em physics lists
	G4PhysListFactory factory;
	g4HadronicList = factory.AvailablePhysLists();
	// em comes with "_", stripping it
	vector<G4String> allEM = factory.AvailablePhysListsEM();
	g4EMList.push_back("STD");
	for (unsigned i = 0; i < allEM.size(); i++) {
		vector<string> emstripped = getStringVectorFromStringWithDelimiter(allEM[i], "_");
		string stripped;

		if (emstripped.size() == 2) stripped = trimSpacesFromString(emstripped[1]);
		else if (emstripped.size() == 1) stripped = emstripped[0];
		else
			continue;

		g4EMList.push_back(stripped);

	}

	G4LossTableManager::Instance();
	defaultCutValue = gemcOpt.optMap["PRODUCTIONCUT"].arg;

	if (fastmcMode > 0) defaultCutValue = 5000;

	cutForGamma = defaultCutValue;
	cutForElectron = defaultCutValue;
	cutForPositron = defaultCutValue;
	cutForProton = defaultCutValue;

	physIngredients = getStringVectorFromStringWithDelimiter(ingredientsList, "+");

	g4EMPhysics = NULL;
	g4ParticleList = NULL;
	g4HadronicPhysics.clear();

	// validateIngredients will also set hadronicPhys, EMPhys, opticalPhys
	if (!validateIngredients()) {
		cout << "  !!! Error: physics ingredients list not valid: >" << ingredientsList << "<" << endl;
		list();

		cout << " Exiting." << endl;
		exit(0);
	} else
		cookPhysics();
}

PhysicsList::~PhysicsList() {

}

void PhysicsList::list() {
	cout << "     > Available hadronic physics list: " << endl;
	for (unsigned i = 0; i < g4HadronicList.size(); i++)
		cout << "       - " << g4HadronicList[i] << endl;
	cout << endl;

	cout << "     > Available EM physics list: " << endl;
	for (unsigned i = 0; i < g4EMList.size(); i++)
		cout << "       - " << g4EMList[i] << endl;

	cout << "   > Optica: optical" << endl;
}

// check the ingredients consistency
bool PhysicsList::validateIngredients() {
	unsigned isHadronicLegit = 0;
	unsigned isEMLegit = 0;
	unsigned isOpticalLegit = 0;
	unsigned isHPLegit = 0;

	G4PhysListFactory factory;
	for (unsigned i = 0; i < physIngredients.size(); i++) {
		string ingredient = trimSpacesFromString(physIngredients[i]);

		if (factory.IsReferencePhysList(ingredient)) {
			isHadronicLegit = 1;
			hadronicPhys = ingredient;
		}

		for (unsigned i = 0; i < g4EMList.size(); i++)
			if (ingredient == g4EMList[i]) {
				isEMLegit = 1;
				EMPhys = ingredient;
			}

		if (ingredient == "Optical") {
			isOpticalLegit = 1;
			opticalPhys = "yes";
		}
		if (ingredient == "HP") {
			HPPhys = "yes";
			isHPLegit = 1;
		}
	}

	if (verbosity > 0) {
		cout << "  >> Physics: " << ingredientsList << endl;
		cout << "   > Hadronic: " << hadronicPhys << endl;
		cout << "   > EM: " << EMPhys << endl;
		cout << "   > HP: " << HPPhys << endl;
		cout << "   > Optical: " << opticalPhys << endl << endl;
	}

	if (physIngredients.size() == (isHadronicLegit + isEMLegit + isOpticalLegit + isHPLegit)) return TRUE;

	return FALSE;
}

void PhysicsList::SetCuts() {

	if (verbosity > 0) {
		cout << "PhysicsList::SetCuts:";
		cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << endl;
	}

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");
	SetCutValue(cutForProton, "proton");

	if (verbosity > 0) DumpCutValuesTable();
}

void PhysicsList::SetCutForGamma(double cut) {
	cutForGamma = cut;
	SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void PhysicsList::SetCutForElectron(double cut) {
	cutForElectron = cut;
	SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void PhysicsList::SetCutForPositron(double cut) {
	cutForPositron = cut;
	SetParticleCuts(cutForPositron, G4Positron::Positron());
}

void PhysicsList::SetCutForProton(double cut) {
	cutForProton = cut;
	SetParticleCuts(cutForProton, G4Proton::Proton());
}

#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
#include "FTFP_BERT_TRV.hh"
#include "FTFP_INCLXX.hh"
#include "FTFP_INCLXX_HP.hh"
#include "FTF_BIC.hh"
#include "LBE.hh"
#include "QBBC.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC_AllHP.hh"
#include "QGSP_FTFP_BERT.hh"
#include "QGS_BIC.hh"
#include "QGSP_INCLXX.hh"
#include "QGSP_INCLXX_HP.hh"
#include "Shielding.hh"
#include "NuBeam.hh"

#include "G4OpticalPhysics.hh"
#include "G4SynchrotronRadiation.hh"
#include "G4SynchrotronRadiationInMat.hh"

#include "G4StepLimiter.hh"

void PhysicsList::cookPhysics() {
	// Particles
	g4ParticleList = new G4DecayPhysics("decays");

	// EM Physics
	if (g4EMPhysics) delete g4EMPhysics;
	if (EMPhys == "STD") g4EMPhysics = new G4EmStandardPhysics();
	else if (EMPhys == "EMV") g4EMPhysics = new G4EmStandardPhysics_option1();
	else if (EMPhys == "EMX") g4EMPhysics = new G4EmStandardPhysics_option2();
	else if (EMPhys == "EMY") g4EMPhysics = new G4EmStandardPhysics_option3();
	else if (EMPhys == "EMZ") g4EMPhysics = new G4EmStandardPhysics_option4();
	else if (EMPhys == "LIV") g4EMPhysics = new G4EmLivermorePhysics();
	else if (EMPhys == "LVP") g4EMPhysics = new G4EmLivermorePolarizedPhysics();
	else if (EMPhys == "PEN") g4EMPhysics = new G4EmPenelopePhysics();
	else if (EMPhys == "LEM") g4EMPhysics = new G4EmLowEPPhysics();
	else if (EMPhys != "none") {
		cout << " !! Wrong EMPhys " << EMPhys << endl << "Exiting." << endl;
		exit(0);
	}

	// Hadronic Physics
	// The lists in FTFP_BERT() functions also include the EM, hadron elastic and so on. Here we separate them.
	// See for example FTFP_BERT.icc
	// This is a general version of
	// Hadr01 example

	// em extra physics always there
	if (hadronicPhys != "none") {
		g4HadronicPhysics.push_back(new G4EmExtraPhysics(verbosity));

		if (HPPhys == "yes") {
			cout << "   >  Loading High Precision Cross Sections... this may take a while..." << endl;
			g4HadronicPhysics.push_back(new G4HadronElasticPhysicsHP(verbosity));
		} else
			g4HadronicPhysics.push_back(new G4HadronElasticPhysics(verbosity));

		// binary cascade, bertini models, or standard
		// ion physics
		if (hadronicPhys.find("BIC") != string::npos) g4HadronicPhysics.push_back(new G4IonBinaryCascadePhysics(verbosity));
		else if (hadronicPhys.find("BERT") != string::npos) g4HadronicPhysics.push_back(new G4IonPhysics(verbosity));

		g4HadronicPhysics.push_back(new G4NeutronTrackingCut(verbosity));
	}

	// adding the hadronic physics list
	// the complete list is in G4PhysListFactory.cc
	// FTF
	if (hadronicPhys == "FTFP_BERT") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsFTFP_BERT());
	} else if (hadronicPhys == "FTFP_BERT_HP") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsFTFP_BERT_HP());
	} else if (hadronicPhys == "FTFP_BERT_TRV") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsFTFP_BERT_TRV());
	} else if (hadronicPhys == "FTF_BIC") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsFTF_BIC());
	}
	// QGSP
	else if (hadronicPhys == "QGSP_BERT") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsQGSP_BERT());
	} else if (hadronicPhys == "QGSP_BERT_HP") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsQGSP_BERT_HP());
	} else if (hadronicPhys == "QGSP_BIC") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsQGSP_BIC());
	} else if (hadronicPhys == "QGSP_BIC_HP") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsQGSP_BIC_HP());
	} else if (hadronicPhys == "QGSP_FTFP_BERT") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsQGSP_FTFP_BERT());
	}
	// Shielding
	else if (hadronicPhys == "Shielding") {
		g4HadronicPhysics.push_back(new G4HadronPhysicsShielding());
	}
	// Others
	else if (hadronicPhys == "none") {
		;
	} else {
		cout << "Wrong hadronicPhys " << hadronicPhys << endl << "Exiting." << endl;
		exit(0);
	}

	// optical physics
	// taken from example: optical/LXe
	if (opticalPhys == "yes") {
		// verbosity is set to zero at the constructor level by default
		// see G4OpticalPhysics.hh
		G4OpticalPhysics *opticalPhysics = new G4OpticalPhysics();
		opticalPhysics->SetWLSTimeProfile("delta");

		g4HadronicPhysics.push_back(opticalPhysics);
	}

}

void PhysicsList::ConstructParticle() {
	muonRadDecay = 0;
	muonRadDecay = gemcOpt.optMap["FORCE_MUON_RADIATIVE_DECAY"].arg;

	g4ParticleList->ConstructParticle();

	if (muonRadDecay) {
		G4DecayTable *MuonPlusDecayTable = new G4DecayTable();
		G4DecayTable *MuonMinusDecayTable = new G4DecayTable();
		MuonPlusDecayTable->Insert(new G4MuonRadiativeDecayChannelWithSpin("mu+", 1.00));
		MuonMinusDecayTable->Insert(new G4MuonRadiativeDecayChannelWithSpin("mu-", 1.00));
		G4MuonPlus::MuonPlusDefinition()->SetDecayTable(MuonPlusDecayTable);
		G4MuonMinus::MuonMinusDefinition()->SetDecayTable(MuonMinusDecayTable);
	}

	/*Here construct the particles for the dark photon model, if requested*/
	if (gemcOpt.optMap["DARK_PHOTON"].args == "no" || gemcOpt.optMap["DARK_MATTER"].args == "no" || gemcOpt.optMap["DARK_COUPLINGS"].args == "no") {
		return;
	}
	vector<string> valuesDP;
	valuesDP = get_info(gemcOpt.optMap["DARK_PHOTON"].args);
	if (valuesDP.size() != 3) {
		cout << " ERROR, DARK_PHOTON should follow with three numbers: mass, 2*J,P quantum numbers, as in 100*MeV,2,-1" << endl;
		exit(1);
	}
	vector<string> valuesDM;
	valuesDM = get_info(gemcOpt.optMap["DARK_MATTER"].args);
	if (valuesDM.size() != 2) {
		cout << " ERROR, DARK_MATTER should follow with two numbers: mass 2*J,  quantum numbers, as in 100*MeV,1" << endl;
		exit(1);
	}
	vector<string> valuesC;
	valuesC = get_info(gemcOpt.optMap["DARK_COUPLINGS"].args);
	if (valuesC.size() != 2) {
		cout << " ERROR, DARK_COUPLINGS should follow with two numbers: eps and alphaD, as in 1E-4,0.1" << endl;
		exit(1);
	}

	G4double darkPhotonMass = get_number(valuesDP[0]);
	G4double darkMatterMass = get_number(valuesDM[0]);
	G4double eps = get_number(valuesC[0]);
	G4double alphaD = get_number(valuesC[1]);

	G4int darkPhotonJ = atoi(valuesDP[1].c_str()); //A.C. no check on the format.
	G4int darkPhotonP = atoi(valuesDP[2].c_str());
	G4int darkMatterJ = atoi(valuesDM[1].c_str()); //A.C. no check on the format.

	if (darkPhotonJ < 0) {
		cout << " ERROR, darkPhoton spin should be >0, you entered: " << darkPhotonJ << endl;
		exit(1);
	}
	if (abs(darkPhotonJ) % 2 != 0) {
		cout << " ERROR, darkPhoton spin should be integer, you entered: " << darkPhotonJ << endl;
		exit(1);
	}

	if (darkMatterJ < 0) {
		cout << " ERROR, darkMatter spin should be >0, you entered: " << darkMatterJ << endl;
		exit(1);
	}

	//First, construct the DM
	G4DarkMatter::DarkMatterDefinition(darkMatterMass, darkMatterJ);
	//Then, construct the DP
	G4DarkPhoton::DarkPhotonDefinition(darkPhotonMass, darkPhotonJ, darkPhotonP, eps, alphaD);

}

void PhysicsList::ConstructProcess() {
	AddTransportation();
	int fastmcMode = gemcOpt.optMap["FASTMCMODE"].arg;
	int synrad = gemcOpt.optMap["SYNRAD"].arg;

	if (fastmcMode % 10 < 2) {
		G4ProcessTable *processTable = G4ProcessTable::GetProcessTable();
		G4VProcess *decay;

		if (g4EMPhysics) g4EMPhysics->ConstructProcess();

		g4ParticleList->ConstructProcess();


		for (size_t i = 0; i < g4HadronicPhysics.size(); i++)

		cout<<"Hadron processes: "<<endl;
		for(size_t i=0; i<g4HadronicPhysics.size(); i++){
		  cout<<i<<": "<<g4HadronicPhysics[i]->GetPhysicsName()<<endl;

			g4HadronicPhysics[i]->ConstructProcess();
		}
 
		
		// sync radiation
		G4SynchrotronRadiation *fSync = nullptr;
		G4SynchrotronRadiationInMat *fSyncMat = nullptr;

		if (synrad == 1) fSync = new G4SynchrotronRadiation();
		else if (synrad == 2) fSyncMat = new G4SynchrotronRadiationInMat();
		//G4AutoDelete::Register(fSync);

		auto theParticleIterator = GetParticleIterator();

		// PhysicsList contains theParticleIterator
		theParticleIterator->reset();

		while ((*theParticleIterator)()) {

			G4ParticleDefinition *particle = theParticleIterator->value();
			G4ProcessManager *pmanager = particle->GetProcessManager();
			string pname = particle->GetParticleName();

			// Adding Step Limiter
			if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0) && (pname != "chargedgeantino")) {
				if (verbosity > 2) cout << "   >  Adding Step Limiter for " << pname << endl;

				pmanager->AddProcess(new G4StepLimiter, -1, -1, 3);
			}

			// G4SynchrotronRadiation if requested
			if (synrad == 1) {

				if (pname == "e-") {
					//electron
					pmanager->AddProcess(fSync, -1, -1, 4);
					pmanager->AddProcess(new G4StepLimiter, -1, -1, 5);

				} else if (pname == "e+") {
					//positron
					pmanager->AddProcess(fSync, -1, -1, 5);
					pmanager->AddProcess(new G4StepLimiter, -1, -1, 6);
				}
			} else if (synrad == 2) {

				if (pname == "e-") {
					//electron
					pmanager->AddProcess(fSyncMat, -1, -1, 4);
					pmanager->AddProcess(new G4StepLimiter, -1, -1, 5);

				} else if (pname == "e+") {
					//positron
					pmanager->AddProcess(fSyncMat, -1, -1, 5);
					pmanager->AddProcess(new G4StepLimiter, -1, -1, 6);
				}
			}

			if (muonRadDecay) {
				theDecayProcess = new G4DecayWithSpin();
				decay = processTable->FindProcess("Decay", particle);
				if (theDecayProcess->IsApplicable(*particle)) {
					if (decay) pmanager->RemoveProcess(decay);
					pmanager->AddProcess(theDecayProcess);
					pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
					pmanager->SetProcessOrderingToLast(theDecayProcess, idxAtRest);
				}
			}
		}


		const G4ParticleDefinition *particle = G4Gamma::Gamma();
		G4ProcessManager *pmanager = particle->GetProcessManager();
		G4GammaConversionToMuons* GCTM = new G4GammaConversionToMuons;
		GCTM -> SetCrossSecFactor(1);
		pmanager->AddDiscreteProcess(GCTM);
		
		
		auto processes=G4Gamma::Gamma()->GetProcessManager()->GetProcessList();
		G4VProcess *proc;
		for (int ii=0;ii<processes->size();ii++){
		  proc=(*processes)[ii];
		  if (proc->GetProcessName()=="photonNuclear"){
		    auto hadr_proc=dynamic_cast<G4HadronInelasticProcess*>(proc);
		    //cout<<"GOT IT: "<<hadr_proc<<" "<<endl;
	      //hadr_proc->BiasCrossSectionByFactor(1E5);
		  }
		}
		
		
		if (gemcOpt.optMap["DARK_PHOTON"].args == "no" || gemcOpt.optMap["DARK_MATTER"].args == "no" || gemcOpt.optMap["DARK_COUPLINGS"].args == "no") {
		  return;
		} else {
		  //A' production
		  vector<string> valuesC;
		  valuesC = get_info(gemcOpt.optMap["DARK_COUPLINGS"].args);
		  if (valuesC.size() != 2) {
		    cout << " ERROR, DARK_COUPLINGS should follow with two numbers: eps and alphaD, as in 1E-4,0.1" << endl;
		    exit(1);
		  }
		  
		  G4double eps = get_number(valuesC[0]);
		  G4double alphaD = get_number(valuesC[1]);
		  G4double eCut = get_number(valuesC[2]);
		  particle = G4Positron::Positron();
		  pmanager = particle->GetProcessManager();
		  G4DarkPhotonAnnihilationProduction* myDarkPhotonAnnihilationProduction = new G4DarkPhotonAnnihilationProduction;
		  myDarkPhotonAnnihilationProduction->SetEps(eps);
		  myDarkPhotonAnnihilationProduction->SetAlphaD(alphaD);
		  myDarkPhotonAnnihilationProduction->SetEnergyCut(eCut);
		  pmanager->AddDiscreteProcess(myDarkPhotonAnnihilationProduction);
		}
	}
}
