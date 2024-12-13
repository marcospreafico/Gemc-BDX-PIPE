// G4 headers
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Trajectory.hh"

// gemc headers
#include "MEventAction.h"
#include "string_utilities.h"
#include "Hit.h"

// mlibrary
#include "frequencySyncSignal.h"

#include <iostream>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/engineIDulong.h"
using namespace CLHEP;

//ROOT headers
#include "TFile.h"
#include "TTree.h"

//EVIO
#include "evioBankUtil.h"
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"

// return original track id of a vector of tid
vector<int> MEventAction::vector_otids(vector<int> tids) {
    vector<int> otids;
    for (unsigned int t = 0; t < tids.size(); t++) {
        if (hierarchy.find(tids[t]) != hierarchy.end()) otids.push_back(hierarchy[tids[t]]);
        else
            otids.push_back(0);
    }
    return otids;
}

// the following functions return the pid, mpid and charges vector<int> from a map<int, TInfos>
vector<int> vector_zint(int size) {
    vector<int> zints;
    for (int t = 0; t < size; t++)
        zints.push_back(0);

    return zints;
}

vector<G4ThreeVector> vector_zthre(int size) {
    vector<G4ThreeVector> zthre;
    for (int t = 0; t < size; t++)
        zthre.push_back(G4ThreeVector(0, 0, 0));

    return zthre;
}

vector<int> vector_mtids(map<int, TInfos> tinfos, vector<int> tids) {
    vector<int> mtids;
    for (unsigned int t = 0; t < tids.size(); t++)
        mtids.push_back(tinfos[tids[t]].mtid);

    return mtids;
}

vector<int> vector_mpids(map<int, TInfos> tinfos, vector<int> tids) {
    vector<int> mpids;
    for (unsigned int t = 0; t < tids.size(); t++)
        mpids.push_back(tinfos[tids[t]].mpid);

    return mpids;
}

vector<G4ThreeVector> vector_mvert(map<int, TInfos> tinfos, vector<int> tids) {
    vector<G4ThreeVector> mvert;
    for (unsigned int t = 0; t < tids.size(); t++)
        mvert.push_back(tinfos[tids[t]].mv);

    return mvert;
}

MEventAction::MEventAction(goptions opts, map<string, double> gpars) {

    hitProcessMap = 0;
    gen_action = 0;
    outputFactoryMap = 0;
    outContainer = 0;
    banksMap = 0;

    gemcOpt = opts;
    hd_msg = gemcOpt.optMap["LOG_MSG"].args + " Event Action: >> ";
    Modulo = (int) gemcOpt.optMap["PRINT_EVENT"].arg;
    VERB = gemcOpt.optMap["BANK_VERBOSITY"].arg;
    catch_v = gemcOpt.optMap["CATCH"].args;
    SAVE_ALL_MOTHERS = (int) gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg;
    gPars = gpars;
    MAXP = (int) gemcOpt.optMap["NGENP"].arg;
    FILTER_HITS = (int) gemcOpt.optMap["FILTER_HITS"].arg;
    rw = runWeights(opts);

    WRITE_ALLRAW = replaceCharInStringWithChars(gemcOpt.optMap["ALLRAWS"].args, ",", "  ");
    WRITE_INTRAW = replaceCharInStringWithChars(gemcOpt.optMap["INTEGRATEDRAW"].args, ",", "  ");
    WRITE_INTDGT = replaceCharInStringWithChars(gemcOpt.optMap["INTEGRATEDDGT"].args, ",", "  ");
    SIGNALVT = replaceCharInStringWithChars(gemcOpt.optMap["SIGNALVT"].args, ",", "  ");
    RFSETUP = replaceCharInStringWithChars(gemcOpt.optMap["RFSETUP"].args, ",", "  ");
    fastMCMode = gemcOpt.optMap["FASTMCMODE"].arg;  // fast mc = 2 will increase prodThreshold and maxStep to 5m

    do_RETRIEVE_RANDOM = false;
    if (((string) gemcOpt.optMap["RETRIEVE_RANDOM"].args) != "no") do_RETRIEVE_RANDOM = true;

    // fastMC mode willset SAVE_ALL_MOTHERS to 1
    // a bit cluncky for now
    if (fastMCMode > 0) SAVE_ALL_MOTHERS = 1;

    tsampling = get_number(get_info(gemcOpt.optMap["TSAMPLING"].args).front());
    nsamplings = get_number(get_info(gemcOpt.optMap["TSAMPLING"].args).back());

    if (SAVE_ALL_MOTHERS > 1) {
        lundOutput = new ofstream("background.dat");
        cout << " > Opening background.dat file to save background particles in LUND format." << endl;
    }

    evtN = gemcOpt.optMap["EVTN"].arg;

    //These are for JPOS trigger
    do_JPOS_TRG = false;
    Eprompt = 0;
    SDprompt = "";
    Eprompt_MIN = 0;
    Eprompt_MAX = 100 * TeV;
    Tprompt_MIN = 0;
    Tprompt_MAX = 100 * ms;

    if (((string) gemcOpt.optMap["JPOS_TRG"].args) != "no") {
        vector<string> cvalues = get_info(gemcOpt.optMap["JPOS_TRG"].args, string(",\""));
        if (cvalues.size() != 5) {
            cout << "ERR: JPOS_TRG must be 5 numbers, separated with comma: SD_name,Eprompt_MIN(with_unit),Eprompt_MAX(with_unit),Tprompt_min(with_unit),Tprompt_max(with_unit)";
            cout << "endl";
            exit(1);
        }
        do_JPOS_TRG = true;
        SDprompt = cvalues[0];
        Eprompt_MIN = get_number(cvalues[1]);
        Eprompt_MAX = get_number(cvalues[2]);
        Tprompt_MIN = get_number(cvalues[3]);
        Tprompt_MAX = get_number(cvalues[4]);

        if (VERB > 3) {
            cout << "MEventAction JPOS_TRG settings:" << endl;
            cout << "SDprompt = " << SDprompt << endl;
            cout << "Eprompt_MIN = " << Eprompt_MIN / GeV << " GeV " << endl;
            cout << "Eprompt_MAX = " << Eprompt_MAX / GeV << " GeV " << endl;
            cout << "Tprompt_MIN = " << Tprompt_MIN / ns << " ns " << endl;
            cout << "Tprompt_MAX = " << Tprompt_MAX / ns << " ns " << endl;

        }
    }

    //These are for JPOS_2 trigger
    do_JPOS_TRG_2 = false;
    Eprompt_2 = 0;
    SDprompt_2 = "";
    Eprompt_MIN_2 = 0;
    Eprompt_MAX_2 = 100 * TeV;
    Tprompt_MIN_2 = 0;
    Tprompt_MAX_2 = 100 * ms;

    if (((string) gemcOpt.optMap["JPOS_TRG_2"].args) != "no") {
        vector<string> cvalues = get_info(gemcOpt.optMap["JPOS_TRG_2"].args, string(",\""));
        if (cvalues.size() != 5) {
            cout << "ERR: JPOS_TRG_2 must be 5 numbers, separated with comma: SD_name,Eprompt_MIN(with_unit),Eprompt_MAX(with_unit),Tprompt_min(with_unit),Tprompt_max(with_unit)";
            cout << "endl";
            exit(1);
        }
        do_JPOS_TRG_2 = true;
        SDprompt_2 = cvalues[0];
        Eprompt_MIN_2 = get_number(cvalues[1]);
        Eprompt_MAX_2 = get_number(cvalues[2]);
        Tprompt_MIN_2 = get_number(cvalues[3]);
        Tprompt_MAX_2 = get_number(cvalues[4]);

        if (VERB > 3) {
            cout << "MEventAction JPOS_TRG_2 settings:" << endl;
            cout << "SDprompt = " << SDprompt_2 << endl;
            cout << "Eprompt_MIN = " << Eprompt_MIN_2 / GeV << " GeV " << endl;
            cout << "Eprompt_MAX = " << Eprompt_MAX_2 / GeV << " GeV " << endl;
            cout << "Tprompt_MIN = " << Tprompt_MIN_2 / ns << " ns " << endl;
            cout << "Tprompt_MAX = " << Tprompt_MAX_2 / ns << " ns " << endl;

        }
    }

    //These are for JPOS_3 trigger (HCAL)
    do_JPOS_TRG_3 = false;
    SDprompt_3 = "";
    Eprompt_MIN_3 = 0;
    Eprompt_MAX_3 = 100 * TeV;
    Tprompt_MIN_3 = 0;
    Tprompt_MAX_3 = 100 * ms;
    for (int ii = 0; ii < 6; ii++)
        Nprompt_bars_3_thr[ii] = 9999;

    if (((string) gemcOpt.optMap["JPOS_TRG_3"].args) != "no") {
        vector<string> cvalues = get_info(gemcOpt.optMap["JPOS_TRG_3"].args, string(",\""));
        if (cvalues.size() != 11) {
            cout << "ERR: JPOS_TRG_3 must be 11 numbers, separated with comma: SD_name,Eprompt_MIN(with_unit),Eprompt_MAX(with_unit),Tprompt_min(with_unit),Tprompt_max(with_unit),THR_SECT_1,THR_SECT_2,THR_SECT_3,THR_SECT_4,THR_SECT_5,THR_SECT_6";
            cout << "endl";
            exit(1);
        }
        do_JPOS_TRG_3 = true;
        SDprompt_3 = cvalues[0];
        Eprompt_MIN_3 = get_number(cvalues[1]);
        Eprompt_MAX_3 = get_number(cvalues[2]);
        Tprompt_MIN_3 = get_number(cvalues[3]);
        Tprompt_MAX_3 = get_number(cvalues[4]);
        for (int ii = 0; ii < 6; ii++)
            Nprompt_bars_3_thr[ii] = get_number(cvalues[5 + ii]);

        if (VERB > 3) {
            cout << "MEventAction JPOS_TRG_3 settings:" << endl;
            cout << "SDprompt = " << SDprompt_3 << endl;
            cout << "Eprompt_MIN = " << Eprompt_MIN_3 / GeV << " GeV " << endl;
            cout << "Eprompt_MAX = " << Eprompt_MAX_3 / GeV << " GeV " << endl;
            cout << "Tprompt_MIN = " << Tprompt_MIN_3 / ns << " ns " << endl;
            cout << "Tprompt_MAX = " << Tprompt_MAX_3 / ns << " ns " << endl;
            cout << "Sector number of bars thr:" << endl;
            for (int ii = 0; ii < 6; ii++)
                cout << ii << " " << Nprompt_bars_3_thr[ii] << endl;
        }
    }

    do_SAVE_RANDOM = false;
    if (gemcOpt.optMap["SAVE_RANDOM"].arg != 0) do_SAVE_RANDOM = true;

}

MEventAction::~MEventAction() {
    if (SAVE_ALL_MOTHERS > 1) lundOutput->close();
}

void MEventAction::BeginOfEventAction(const G4Event *evt) {

    if (gen_action->isFileOpen() == false) {
        G4RunManager *runManager = G4RunManager::GetRunManager();
        ;
        runManager->AbortRun();
        cout << " No more events in the input file." << endl;
        return;
    }

    rw.getRunNumber(evtN);
    bgMap.clear();

    //A.C. this part was moved to the MPrimaryGeneratorAction, since random numbers are already used there, and that function is called before
    /* if (do_RETRIEVE_RANDOM) {
     retrieveRandom();
     }*/

   // header.clear();
    //A.C. this part was moved to the MPrimaryGeneratorAction, since random numbers are already used there, and that function is called before
  /*  if (do_SAVE_RANDOM) {
      this->saveRandom();
    }*/

    if (evtN % Modulo == 0) {
        cout << hd_msg << " Begin of event " << evtN << "  Run Number: " << rw.runNo;
        if (rw.isNewRun) cout << " (new) ";
        cout << endl;
        cout << hd_msg << " Random Number: " << G4UniformRand() << endl;
        // CLHEP::HepRandom::showEngineStatus();

    }

    //JPOS_CRS part
    Eprompt = 0;
    Eprompt_2 = 0;
    for (int ii = 0; ii < 6; ii++) {
        for (int jj = 0; jj < kMaxSector; jj++) {
            for (int kk = 0; kk < kMaxChannel; kk++) {
                Eprompt_bars_3[ii][jj][kk] = 0;
            }
        }
    }

}

void MEventAction::EndOfEventAction(const G4Event *evt) {
    if (gen_action->isFileOpen() == false) {
        return;
    }

    MHitCollection *MHC;
    int nhits;

    // if FILTER_HITS is set, checking if there are any hits
    if (FILTER_HITS) {
        int anyHit = 0;
        for (map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it != SeDe_Map.end(); it++) {
            MHC = it->second->GetMHitCollection();
            if (MHC) anyHit += MHC->GetSize();
        }

        // stop here if there are no hits and FILTER_HITS is set
        if (anyHit == 0) return;
    }

    if (evtN % Modulo == 0) cout << hd_msg << " Starting Event Action Routine " << evtN << "  Run Number: " << rw.runNo << endl;

    if (do_JPOS_TRG || do_JPOS_TRG_2 || do_JPOS_TRG_3) {
        bool saveEvent = true;
        if (VERB > 3) {
            cout << "EndOfEventAction " << evtN << endl;
            cout << "Eprompt JPOS_TRG: " << Eprompt / GeV << " GeV " << endl;
            cout << "Eprompt JPOS_TRG_2: " << Eprompt_2 / GeV << " GeV " << endl;
        }
        if (do_JPOS_TRG) {
            if ((Eprompt < Eprompt_MIN) || (Eprompt > Eprompt_MAX)) {
                saveEvent = false;
            }
        }
        if (do_JPOS_TRG_2) {
            if ((Eprompt_2 < Eprompt_MIN_2) || (Eprompt_2 > Eprompt_MAX_2)) {
                saveEvent = false;
            }
        }

        if (do_JPOS_TRG_3) {
            bool tmpFlag = true;
            for (int ii = 0; ii < 6; ii++) {
                Nprompt_bars_3[ii] = 0;
                for (int jj = 0; jj < kMaxSector; jj++) {
                    for (int kk = 0; kk < kMaxChannel; kk++) {
                        if ((Eprompt_bars_3[ii][jj][kk] >= Eprompt_MIN_3) && (Eprompt_bars_3[ii][jj][kk] <= Eprompt_MAX_3)) {
                            Nprompt_bars_3[ii] += 1;
                            if (VERB > 3) cout << "JPOS_TRG_3 sect: " << ii + 1 << " lay: " << jj + 1 << " ch: " << kk + 1 << " " << Eprompt_bars_3[ii][jj][kk] << endl;
                        }
                    }
                }
                if (Nprompt_bars_3[ii] > Nprompt_bars_3_thr[ii]) tmpFlag = false;
            }
            if (tmpFlag == false) { //if at least one of the 6 sectors has a number of bars ON greater than the THR do NOT save
                saveEvent = false;
            }
        }
        if (!saveEvent) {
            if (VERB > 3) {
                cout << " DO NOT SAVE" << endl;
            }
            evtN++;
            return;
        }
    }

    // building the tracks set database with all the tracks in all the hits
    // if SAVE_ALL_MOTHERS is set
    set<int> track_db;
    if (SAVE_ALL_MOTHERS) for (map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it != SeDe_Map.end(); it++) {
        string hitType;

        MHC = it->second->GetMHitCollection();
        nhits = MHC->GetSize();

        for (int h = 0; h < nhits; h++) {
            vector<int> tids = (*MHC)[h]->GetTIds();

            for (unsigned int t = 0; t < tids.size(); t++) {
                track_db.insert(tids[t]);
            }
            if (SAVE_ALL_MOTHERS > 1) {
                vector<int> pids = (*MHC)[h]->GetPIDs();
                // getting the position of the hit, not vertex of track
        vector<G4ThreeVector> vtxs = (*MHC)[h]->GetPos();
        vector<G4ThreeVector> mmts = (*MHC)[h]->GetMoms();
        vector<double> tims = (*MHC)[h]->GetTime();
        // only put the first step of a particular track
        // (don't fill if track exist already)
        for (unsigned int t = 0; t < tids.size(); t++) {
            if (bgMap.find(tids[t]) == bgMap.end()) bgMap[tids[t]] = BGParts(pids[t], tims[t], vtxs[t], mmts[t]);
        }
    }
}
}

    // now filling the map of tinfos with tracks infos from the track_db database
    // this won't get the mother particle infos except for their track ID
    map<int, TInfos> tinfos;
    set<int>::iterator it;

    // the container is full only if /tracking/storeTrajectory 2
    G4TrajectoryContainer *trajectoryContainer;

    if (SAVE_ALL_MOTHERS) {
        trajectoryContainer = evt->GetTrajectoryContainer();
        momDaughter.clear();

        if (VERB > 3) cout << " >> Total number of tracks " << trajectoryContainer->size() << endl;

        while (trajectoryContainer && track_db.size()) {
            // looping over all tracks
            for (unsigned int i = 0; i < trajectoryContainer->size(); i++) {
                // cout << " index track " << i << endl;
                G4Trajectory *trj = (G4Trajectory*) (*(evt->GetTrajectoryContainer()))[i];
                int tid = trj->GetTrackID();

                // adding track in mom daughter relationship
                // for all tracks
                if (momDaughter.find(tid) == momDaughter.end()) momDaughter[tid] = trj->GetParentID();

                // is this track involved in a hit?
                // if yes, add it to the track map db
                it = track_db.find(tid);
                if (it != track_db.end()) {
                    int mtid = trj->GetParentID();
                    tinfos[tid] = TInfos(mtid);
                    // cout << " At map level: " << tid << " " <<mtid << " " << pid << "   charge: " << q << endl;
                    // remove the track from the database so we don't have to loop over it next time
                    track_db.erase(it);
                }
            }
        }

        // building the hierarchy map
        // for all secondary tracks
        for (map<int, int>::iterator itm = momDaughter.begin(); itm != momDaughter.end(); itm++) {
            int ancestor = itm->first;
            if (momDaughter[ancestor] == 0) hierarchy[itm->first] = itm->first;

            while (momDaughter[ancestor] != 0) {
                hierarchy[itm->first] = momDaughter[ancestor];
                ancestor = momDaughter[ancestor];
            }
        }

        // now accessing the mother particle infos
        for (map<int, TInfos>::iterator itm = tinfos.begin(); itm != tinfos.end(); itm++) {
            int mtid = (*itm).second.mtid;
            // looking for mtid infos in the trajectoryContainer
            for (unsigned int i = 0; i < trajectoryContainer->size() && mtid != 0; i++) {
                G4Trajectory *trj = (G4Trajectory*) (*(evt->GetTrajectoryContainer()))[i];
                int tid = trj->GetTrackID();
                if (tid == mtid) {
                    tinfos[(*itm).first].mpid = trj->GetPDGEncoding();
                    tinfos[(*itm).first].mv = trj->GetPoint(0)->GetPosition();
                }
            }
        }
        // removing daughter tracks if mother also gave hits
        if (SAVE_ALL_MOTHERS > 2) {
            vector<int> bgtIDs;
            for (map<int, BGParts>::iterator it = bgMap.begin(); it != bgMap.end(); it++)
                bgtIDs.push_back(it->first);

            for (unsigned i = 0; i < bgtIDs.size(); i++) {
                int daughter = bgtIDs[i];
                int momCheck = 0;
                if (momDaughter.find(daughter) != momDaughter.end()) momCheck = momDaughter[daughter];

                while (momCheck != 0) {
                    // mom has a hit
                    if (bgMap.find(momCheck) != bgMap.end()) {
                        // so if daughter is found, delete it
                        // daughter may already be deleted before
                        if (bgMap.find(daughter) != bgMap.end()) bgMap.erase(bgMap.find(daughter));
                    }
                    // going up one generation
                    if (momDaughter.find(momCheck) != momDaughter.end()) momCheck = momDaughter[momCheck];
                    else
                        momCheck = 0;
                }
            }
        }
    }

    // Making sure the output routine exists in the ProcessOutput map
    // If no output selected, or HitProcess not found, end current event
    map<string, outputFactoryInMap>::iterator ito = outputFactoryMap->find(outContainer->outType);
    if (ito == outputFactoryMap->end()) {
        if (outContainer->outType != "no") cout << hd_msg << " Warning: output type <" << outContainer->outType << "> is not registered in the outContainerput Factory. " << endl << "      This event will not be written out." << endl;
        evtN++;
        return;
    }
    outputFactory *processOutputFactory = getOutputFactory(outputFactoryMap, outContainer->outType);

    // Header Bank contains event number
    // Need to change this to read DB header bank
    header["runNo"] = rw.runNo;
    header["evn"] = evtN;
    header["evn_type"] = -1;  // physics event. Negative is MonteCarlo event
    header["beamPol"] = gen_action->getBeamPol();

    // write event header bank
    processOutputFactory->writeHeader(outContainer, header, getBankFromMap("header", banksMap));

    // write event header bank
    //(A.C. this is what is used by BDX for event weight, and in general to store LUND header!)
    map<string, double> userHeader;
    for (unsigned i = 0; i < gen_action->headerUserDefined.size(); i++) {
        string tmp = "userVar";
        if (i < 9) tmp += "00";
        else if (i < 99) tmp += "0";

        tmp += to_string(i + 1);

        userHeader[tmp] = gen_action->headerUserDefined[i];
    }
    processOutputFactory->writeUserInfoseHeader(outContainer, userHeader);

    // write RF bank if present
    // do not write in FASTMC mode
    if (RFSETUP != "no" && fastMCMode == 0) {
        vector<string> rfvalues = getStringVectorFromString(RFSETUP);

        // getting time window
        string rfsetup = to_string(gen_action->getTimeWindow()) + " ";

        // getting start time of the event
        rfsetup += to_string(gen_action->getStartTime()) + " ";

        for (unsigned i = 0; i < rfvalues.size(); i++)
            rfsetup += rfvalues[i] + " ";

        FrequencySyncSignal rfs(rfsetup);
        processOutputFactory->writeRFSignal(outContainer, rfs, getBankFromMap("rf", banksMap));

        if (VERB > 1) cout << rfs << endl;
    }

    // Getting Generated Particles info
    // Are these loops necessary, revisit later1
    vector<generatedParticle> MPrimaries;
    for (int pv = 0; pv < evt->GetNumberOfPrimaryVertex() && pv < MAXP; pv++) {
        generatedParticle Mparticle;
        G4PrimaryVertex *MPV = evt->GetPrimaryVertex(pv);
        Mparticle.vertex = MPV->GetPosition();
        double thisTime = MPV->GetT0();
        int thisMult = MPV->GetNumberOfParticle();

        for (int pp = 0; pp < MPV->GetNumberOfParticle() && pv < MAXP; pp++) {
            G4PrimaryParticle *PP = MPV->GetPrimary(pp);
            Mparticle.momentum = PP->GetMomentum();
            Mparticle.PID = PP->GetPDGcode();
            Mparticle.time = thisTime;
            Mparticle.multiplicity = thisMult;
        }
        MPrimaries.push_back(Mparticle);
    }

    if (SAVE_ALL_MOTHERS > 1) saveBGPartsToLund();

    for (map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it != SeDe_Map.end(); it++) {
        MHC = it->second->GetMHitCollection();
        nhits = MHC->GetSize();

        // The same ProcessHit Routine must apply to all the hits  in this HitCollection.
        // Instantiating the ProcessHitRoutine only once for the first hit.
        if (nhits) {
            // the bank idtag is the one that corresponds to the hitType
            //MHit* aHit = (*MHC)[0];
            string vname = (*MHC)[0]->GetDetector().name;
            string hitType = it->second->GetDetectorHitType(vname);

            HitProcess *hitProcessRoutine = getHitProcess(hitProcessMap, hitType, vname);
            if (!hitProcessRoutine) return;

            if (fastMCMode == 0 || fastMCMode > 9) hitProcessRoutine->init(hitType, gemcOpt, gPars);

            bool WRITE_TRUE_INTEGRATED = 0;
            bool WRITE_TRUE_ALL = 0;
            if (WRITE_INTRAW.find(hitType) != string::npos) WRITE_TRUE_INTEGRATED = 1;
            if (WRITE_ALLRAW.find(hitType) != string::npos) WRITE_TRUE_ALL = 1;

            vector<hitOutput> allRawOutput;

            // creating summary information for each generated particle
            for (unsigned pi = 0; pi < MPrimaries.size(); pi++) {
                MPrimaries[pi].pSum.push_back(summaryForParticle("na"));
                if (fastMCMode > 0) MPrimaries[pi].fastMC.push_back(fastMCForParticle("na"));
            }

            for (int h = 0; h < nhits; h++) {
                MHit *aHit = (*MHC)[h];
                if (aHit->isElectronicNoise) continue;

                hitOutput thisHitOutput;

                // mother particle infos
                if (SAVE_ALL_MOTHERS) {
                    // setting track infos before processing the hit
                    vector<int> tids = aHit->GetTIds();
                    vector<int> otids = vector_otids(tids);
                    aHit->SetmTrackIds(vector_mtids(tinfos, tids));
                    aHit->SetoTrackIds(otids);
                    aHit->SetmPIDs(vector_mpids(tinfos, tids));
                    aHit->SetmVerts(vector_mvert(tinfos, tids));

                    // for every particle initializing a vector
                    // the index is the primary particle index
                    // the int will increase for each step
                    // if the int > 0
                    // then we count it as ONE hit by the track
                    vector<int> hitByPrimary;
                    for (unsigned pi = 0; pi < MPrimaries.size(); pi++)
                        hitByPrimary.push_back(0);

                    // all these vector have the same length.
                    for (unsigned pi = 0; pi < MPrimaries.size(); pi++) {
                        vector<double> edeps = aHit->GetEdep();
                        vector<double> times = aHit->GetTime();
                        MPrimaries[pi].pSum.back().nphe = aHit->GetTIds().size();
                        if (fastMCMode > 0) {
                            MPrimaries[pi].fastMC.back().pOrig = aHit->GetMom();
                            MPrimaries[pi].fastMC.back().pSmear = hitProcessRoutine->psmear(aHit->GetMom());
                        }
                        for (unsigned ss = 0; ss < edeps.size(); ss++) {
                            if (otids[ss] == (int) pi + 1) {
                                MPrimaries[pi].pSum.back().etot += edeps[ss];
                                hitByPrimary[pi]++;
                                // getting fastest time - should we put threshold here?
                                if (MPrimaries[pi].pSum.back().t < 0 || MPrimaries[pi].pSum.back().t > times[ss]) MPrimaries[pi].pSum.back().t = times[ss];

                            }
                        }

                        if (hitByPrimary[pi]) MPrimaries[pi].pSum.back().stat++;

                        if (MPrimaries[pi].pSum.back().etot > 0 || MPrimaries[pi].pSum.back().nphe > 0) MPrimaries[pi].pSum.back().dname = hitType;
                    }
                } else {
                    // filling mother infos with zeros
                    int thisHitSize = aHit->GetId().size();
                    vector<int> zint = vector_zint(thisHitSize);
                    vector<G4ThreeVector> zthre = vector_zthre(thisHitSize);
                    aHit->SetoTrackIds(zint);
                    aHit->SetmTrackIds(zint);
                    aHit->SetmPIDs(zint);
                    aHit->SetmVerts(zthre);
                }

                if (fastMCMode == 0 || fastMCMode > 9) thisHitOutput.setRaws(hitProcessRoutine->integrateRaw(aHit, h + 1, WRITE_TRUE_INTEGRATED));

                if (WRITE_TRUE_ALL && (fastMCMode == 0 || fastMCMode > 9)) thisHitOutput.setAllRaws(hitProcessRoutine->allRaws(aHit, h + 1));

                allRawOutput.push_back(thisHitOutput);

                string vname = aHit->GetId()[aHit->GetId().size() - 1].name;
                if (VERB > 4 || vname.find(catch_v) != string::npos) {
                    cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
                    cout << aHit->GetId();
                    double Etot = 0;
                    for (unsigned int e = 0; e < aHit->GetPos().size(); e++)
                        Etot = Etot + aHit->GetEdep()[e];
                    cout << "   Total energy deposited: " << Etot / MeV << " MeV" << endl;
                }
            }

            // geant4 integrated raw information
            // by default they are all DISABLED
            // user can enable them one by one
            // using the INTEGRATEDRAW option
            if (WRITE_TRUE_INTEGRATED) processOutputFactory->writeG4RawIntegrated(outContainer, allRawOutput, hitType, banksMap);

            // geant4 all raw information
            // by default they are all DISABLED
            // user can enable them one by one
            // using the ALLRAWS option
            if (WRITE_TRUE_ALL) processOutputFactory->writeG4RawAll(outContainer, allRawOutput, hitType, banksMap);

            if (VERB > 4) for (unsigned pi = 0; pi < MPrimaries.size(); pi++) {
                cout << " Particle " << pi + 1 << " has " << MPrimaries[pi].pSum.size() << " particle summaries:" << endl;
                for (unsigned ss = 0; ss < MPrimaries[pi].pSum.size(); ss++) {
                    cout << " \t det: " << MPrimaries[pi].pSum[ss].dname << "  Etot: " << MPrimaries[pi].pSum[ss].etot << "  time: " << MPrimaries[pi].pSum[ss].t << endl;
                }
                cout << endl;

            }

            // geant4 integrated digitized information
            // by default they are all ENABLED
            // user can disable them one by one
            // using the INTEGRATEDDGT option
            // for FASTMC mode, do not digitize the info
            if (WRITE_INTDGT.find(hitType) == string::npos && (fastMCMode == 0 || fastMCMode > 9)) {
                hitProcessRoutine->initWithRunNumber(rw.runNo);

                vector<hitOutput> allDgtOutput;
                for (int h = 0; h < nhits; h++) {

                    hitOutput thisHitOutput;
                    MHit *aHit = (*MHC)[h];

                    thisHitOutput.setDgtz(hitProcessRoutine->integrateDgt(aHit, h + 1));
                    allDgtOutput.push_back(thisHitOutput);

                    string vname = aHit->GetId()[aHit->GetId().size() - 1].name;
                    if (VERB > 4 || vname.find(catch_v) != string::npos) {
                        cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
                        cout << aHit->GetId();
                        double Etot = 0;
                        for (unsigned int e = 0; e < aHit->GetPos().size(); e++)
                            Etot = Etot + aHit->GetEdep()[e];
                        cout << "   Total energy deposited: " << Etot / MeV << " MeV" << endl;
                    }
                }
                processOutputFactory->writeG4DgtIntegrated(outContainer, allDgtOutput, hitType, banksMap);

            } // end of geant4 integrated digitized information

            // geant4 voltage versus time
            // by default they are all DISABLED
            // user can enable them one by one
            // using the SIGNALVT option
            if (SIGNALVT.find(hitType) != string::npos) {
                vector<hitOutput> allVTOutput;

                for (int h = 0; h < nhits; h++) {

                    hitOutput thisHitOutput;

                    MHit *aHit = (*MHC)[h];

                    // process each step to produce a charge/time digitized information / step
                    thisHitOutput.setChargeTime(hitProcessRoutine->chargeTime(aHit, h));

                    vector<double> stepTimes = thisHitOutput.getChargeTime()[3];
                    vector<double> stepCharges = thisHitOutput.getChargeTime()[2];
                    vector<double> hardware = thisHitOutput.getChargeTime()[5];

                    map<int, int> vSignal;

                    // crate, slot, channels as from translation table
                    vSignal[0] = hardware[0];
                    vSignal[1] = hardware[1];
                    vSignal[2] = hardware[2];

                    for (unsigned ts = 0; ts < nsamplings; ts++) {
                        double forTime = ts * tsampling;
                        double voltage = 0;

                        // create the voltage output based on the hit process
                        // routine voltage(double charge, double time, double forTime)
                        for (unsigned s = 0; s < stepTimes.size(); s++) {

                            double stepTime = stepTimes[s];
                            double stepCharge = stepCharges[s];

                            voltage += hitProcessRoutine->voltage(stepCharge, stepTime, forTime);

                            // cout << " hit " << h <<    "step " << s << "  time: " << stepTime
                            // << "   charge " << stepCharge << "  voltage " << voltage << "  for time bunch " << ts << endl;
                        }
                        // need conversion factor from double to int
                        // the first 3 entries are crate/slot/channels above
                        vSignal[ts + 3] = (int) voltage;
                    }
                    thisHitOutput.createQuantumS(vSignal);

                    allVTOutput.push_back(thisHitOutput);

                    // this is not written out yet

                    string vname = aHit->GetId()[aHit->GetId().size() - 1].name;

                    if (VERB > 4 || vname.find(catch_v) != string::npos) {
                        cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
                        cout << aHit->GetId();
                        double Etot = 0;
                        for (unsigned int e = 0; e < aHit->GetPos().size(); e++)
                            Etot = Etot + aHit->GetEdep()[e];
                        cout << "   Total energy deposited: " << Etot / MeV << " MeV" << endl;
                    }
                }
                processOutputFactory->writeChargeTime(outContainer, allVTOutput, hitType, banksMap);
                processOutputFactory->writeFADCMode1(outContainer, allVTOutput);
            }

            delete hitProcessRoutine;
        }
    }
    // writing out generated particle infos
    processOutputFactory->writeGenerated(outContainer, MPrimaries, banksMap, gen_action->userInfo);

    processOutputFactory->writeEvent(outContainer);
    delete processOutputFactory;

    if (evtN % Modulo == 0) cout << hd_msg << " End of Event " << evtN << " Routine..." << endl << endl;

    // Increase event number. Notice: this is different than evt->GetEventID()
    evtN++;

    return;
}

void MEventAction::saveBGPartsToLund() {
    // for lund format see documentation at gemc.jlab.org:
    // https://gemc.jlab.org/gemc/html/documentation/generator/lund.html

    // this should work also for no hits, map size will be zero.

    *lundOutput << (int) bgMap.size() << "\t" << evtN << "\t 0 0 0 0 0 0 0 0 " << endl;

    int i = 1;
    for (map<int, BGParts>::iterator it = bgMap.begin(); it != bgMap.end(); it++)
        *lundOutput << i++ << "\t0\t1\t" << it->second.pid << "\t0\t" << it->first << "\t" << it->second.p.x() / GeV << "\t" << it->second.p.y() / GeV << "\t" << it->second.p.z() / GeV << "\t" << it->second.time << "\t0\t" << it->second.v.x() / cm << "\t" << it->second.v.y() / cm
                << "\t" << it->second.v.z() / cm << endl;
}

void MEventAction::saveRandom() {
    header.clear();
    HepRandomEngine *rndm = CLHEP::HepRandom::getTheEngine();
    CLHEP::MTwistEngine *rndm_engine_mtwist = NULL;
    rndm_engine_mtwist = dynamic_cast<CLHEP::MTwistEngine*>(rndm);
    if (rndm_engine_mtwist != NULL) {
        std::vector<unsigned long> v = rndm_engine_mtwist->put(); //626 entries: engineIDulong<MTwistEngine>(), then the 624 numbers in array mt, then the count624 variable
        header["userVarUL0000 MTwistEngine"] = 1;
        header["userVarUL0001"] = CLHEP::HepRandom::getTheSeed();
        for (unsigned i = 0; i < v.size(); i++) {
            string tmp = "userVarUL";
            if (i < 8) tmp += "000";
            else if (i < 98) tmp += "00";
            else if (i < 998) tmp += "0";
            tmp += to_string(i + 2);
            header[tmp] = v[i];
        }
    }
}

void MEventAction::retrieveRandom() {
    vector<string> cvalues = get_info(gemcOpt.optMap["RETRIEVE_RANDOM"].args, string(",\""));
    if (cvalues.size() != 3) {
        cout << "ERR: RETRIEVE_RANDOM must be 3 value, separated with comma: file(.root or .evio),runNumber,eventNumber";
        cout << "endl";
        exit(1);
    }

    string file = cvalues[0];
    int runN = atoi(cvalues[1].c_str());
    int eventN = atoi(cvalues[2].c_str());

    bool doRoot = false;
    bool doEvio = false;

    if (file.find(".root") != string::npos) {
        cout << "MPrimaryGeneratorAction::retrieveRandom() ROOT FILE: " << file << endl;
        doRoot = true;
    } else if (file.find(".evio") != string::npos) {
        cout << "MPrimaryGeneratorAction::retrieveRandom() EVIO FILE: " << file << endl;
        doEvio = true;
    } else {
        cout << "MPrimaryGeneratorAction::retrieveRandom() FILE FORMAT NOT RECOGNIZED: " << file << endl;
        exit(1);
    }

    if (doRoot) retrieveRandomFromRoot(file, runN, eventN);
    else if (doEvio) retrieveRandomFromEvio(file, runN, eventN);
}

void MEventAction::retrieveRandomFromRoot(string file, int runN, int eventN) {
    TFile *f = new TFile(file.c_str());
    TTree *t = (TTree*) f->Get("header");
    if (t != NULL) {
        vector<string> *string_tree = 0;
        vector<double> *runN_tree = 0;
        vector<double> *eventN_tree = 0;

        t->SetBranchAddress("runNo", &runN_tree);
        t->SetBranchAddress("evn", &eventN_tree);
        t->SetBranchAddress("user", &string_tree);
        for (int ii = 0; ii < t->GetEntries(); ii++) {
            t->GetEntry(ii);
            if (((int) ((*runN_tree)[0]) == runN) && ((int) ((*eventN_tree)[0]) == eventN)) {
                cout << "MPrimaryGeneratorAction::retrieveRandomFromRoot() found: " << runN << " " << eventN << " at pos: " << ii << endl;
                retreiveRandomFromString((*string_tree)[0]);
                delete f;
                return;
            }
        }
        cout << "MPrimaryGeneratorAction::retrieveRandomFromRoot() runN and eventN not found " << runN << " " << eventN << endl;
        delete f;
        exit(1);
        return;
    }
    cout << "MPrimaryGeneratorAction::retrieveRandomFromRoot() runN and eventN not found " << runN << " " << eventN << endl;
    exit(1);
    delete f;
    return;
}

void MEventAction::retrieveRandomFromEvio(string file, int runN, int eventN) {

    evioFileChannel *chan = new evioFileChannel(file.c_str(), "r", 3000000);
    chan->open();
    int runN_file, eventN_file;
    string data;
    while (chan->read()) {
        evioDOMTree *EDT = new evioDOMTree(chan);
        evio::evioDOMNodeListP fullList = EDT->getNodeList();
        evio::evioDOMNodeList::const_iterator iter;
        for (iter = fullList->begin(); iter != fullList->end(); iter++) {
            if (((*iter)->tag == HEADER_TAG) && ((*iter)->num == HEADER_NUM_RUNNO)) {
                const evio::evioDOMLeafNode<long> *leaf = static_cast<const evio::evioDOMLeafNode<long>*>(*iter);
                runN_file = leaf->data[0];
            } else if (((*iter)->tag == HEADER_TAG) && ((*iter)->num == HEADER_NUM_EVN)) {
                const evio::evioDOMLeafNode<long> *leaf = static_cast<const evio::evioDOMLeafNode<long>*>(*iter);
                eventN_file = leaf->data[0];
            } else if (((*iter)->tag == HEADER_TAG) && ((*iter)->num == HEADER_NUM_DATA)) {
                const evio::evioDOMLeafNode<string> *leaf = static_cast<const evio::evioDOMLeafNode<string>*>(*iter);
                data = leaf->data[0];
            }
        }
        if ((runN_file == runN) && (eventN_file == eventN)) {
            retreiveRandomFromString(data);
            delete EDT;
            delete chan;
            return;
        }
        delete EDT;
    }

    cout << "MPrimaryGeneratorAction::retrieveRandomFromEvio() runN and eventN not found " << runN << " " << eventN << endl;
    exit(1);
    delete chan;

}

void MEventAction::retreiveRandomFromString(string data) {
    cout << "MPrimaryGeneratorAction:retreiveRandomFromString " << endl;

    /*A.C. according to the code in MEventAction.cc
     the string is organized as:
     userVarUL0000 TYPE_OF_RANDOM_NUMBER_GENERATOR 1
     userVarUL0001 data specific to the random_number_generator
     userVarUL0002
     ...
     */
    istringstream strm(data);
    string user, type;
    int tmpOne;
    strm >> user >> type >> tmpOne;

    if (type == "MTwistEngine") {
        cout << "MPrimaryGeneratorAction:retreiveRandomFromString MTwistEngine " << type << " " << engineIDulong<MTwistEngine>() << endl;

        HepRandomEngine *rndm = CLHEP::HepRandom::getTheEngine();
        CLHEP::MTwistEngine *rndm_engine_mtwist = NULL;
        rndm_engine_mtwist = dynamic_cast<CLHEP::MTwistEngine*>(rndm);
        if (rndm_engine_mtwist == NULL) {
            cout << " The current engine is NOT MTwistEngine " << endl;
            return;
        } else {
            unsigned long seed;
            unsigned long udata;
            strm >> user >> seed;
            std::vector<unsigned long> data;
            for (int ii = 0; ii < 626; ii++) {
                strm >> user >> udata;
                data.push_back(udata);
            }
            CLHEP::HepRandom::setTheSeed(seed);
            rndm_engine_mtwist->get(data);
        }
    } else {
        cout << "MPrimaryGeneratorAction:retreiveRandomFromString TYPE NOT RECOGNIZED " << type << endl;
    }
}

