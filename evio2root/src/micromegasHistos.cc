// EVIO headers
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"
using namespace evio;

// gemc headers
#include "options.h"
#include "gbank.h"

// C++ headers
#include <string>
#include <iostream>
#include <vector>
using namespace std;

// banks header
#include "banks.h"

// ROOT Headers
#include "TFile.h"
#include "TTree.h"


int main(int argc, char **argv)
{
	// argument is input file name.
	try
	{
		// loading options
		goptions gemcOpt;
		gemcOpt.setGoptions();
		gemcOpt.setOptMap(argc, argv);
		
		// reading input	  
		string inputfile = gemcOpt.optMap["INPUTF"].args ;
		string bank = gemcOpt.optMap["B"].args ;		
		
		double verbosity = gemcOpt.optMap["BANK_VERBOSITY"].arg;
		
		// loading veriables definitions from MYSQL db
		map<string, string> bank_factories;
// 		bank_factories[bank] = "MYSQL Original";
		bank_factories[bank] = "TEXT original";		
		
		vector<string> whichSystems;
		whichSystems.push_back(bank);
		map<string, gBank> banksMap = read_banks(gemcOpt, bank_factories, whichSystems);
		
		// ROOT file for output; defining tree
		// ROOT does not recognize vector<float> but vector<doubles> are fine
		string outputfile = inputfile.substr(0,inputfile.rfind("."))+".root";
		TFile *f = new TFile(outputfile.c_str(),"RECREATE");

		int evn;
		vector<int> layer;
		vector<int> sector;
		vector<int> strip;		
		vector<double> pid;
		vector<double> mpid;
		vector<double> tid;
		vector<double> mtid;
		vector<double> trackE;
		vector<double> totEdep;
		vector<double> gX;
		vector<double> gY;
		vector<double> gZ;
		vector<double> lX;
		vector<double> lY;
		vector<double> lZ;
		vector<double> px;
		vector<double> py;
		vector<double> pz;
		vector<double> vx;
		vector<double> vy;
		vector<double> vz;
		vector<double> mvx;
		vector<double> mvy;
		vector<double> mvz;
		vector<double> t;
		
		// creating the bank tree
		TTree *bankT = new TTree((bank+"T").c_str(),(bank+" Tree").c_str());
		bankT->Branch("evn" ,   &evn  );
		
		bankT->Branch("layer", &layer);
		bankT->Branch("sector", &sector);
		bankT->Branch("strip", &strip);
		
		bankT->Branch("pid",     &pid);
		bankT->Branch("mpid",    &mpid);
		bankT->Branch("tid",     &tid);
		bankT->Branch("mtid",    &mtid);
		bankT->Branch("trackE",  &trackE);
		bankT->Branch("totEdep", &totEdep);
		
		bankT->Branch("gX" ,     &gX);
		bankT->Branch("gY" ,     &gY);
		bankT->Branch("gZ" ,     &gZ);
		bankT->Branch("lX" ,     &lX);
		bankT->Branch("lY" ,     &lY);
		bankT->Branch("lZ" ,     &lZ);
				
		bankT->Branch("px" ,     &px);
		bankT->Branch("py" ,     &py);
		bankT->Branch("pz" ,     &pz);
		bankT->Branch("vx" ,     &vx);
		bankT->Branch("vy" ,     &vy);
		bankT->Branch("vz" ,     &vz);
		bankT->Branch("mvx" ,    &mvx);
		bankT->Branch("mvy" ,    &mvy);
		bankT->Branch("mvz" ,    &mvz);
		bankT->Branch("t" ,      &t);

		// creating the generated particle info tree

		vector<double> gen_pid;
		vector<double> gen_p;
		vector<double> gen_theta;
		vector<double> gen_phi;
		vector<double> gen_vx;
		vector<double> gen_vy;
		vector<double> gen_vz;

		
		TTree *genT = new TTree("genT", "generated particles Tree");
		genT->Branch("pid"  ,   &gen_pid);
		genT->Branch("p"    ,   &gen_p);
		genT->Branch("theta",   &gen_theta);
		genT->Branch("phi"  ,   &gen_phi);
		genT->Branch("vx"   ,   &gen_vx);
		genT->Branch("vy"   ,   &gen_vy);
		genT->Branch("vz"   ,   &gen_vz);
		

		evioFileChannel *chan = new evioFileChannel(inputfile, "r", 3000000);
		chan->open();
		while(chan->read())
		{
			evioDOMTree EDT(chan);
		
			// read hit banks
			for(map<string, gBank>::iterator it = banksMap.begin(); it != banksMap.end(); it++)
			{
// 				int counter=0;					
				// only considering data banks (bankType.size() > 0) different from "raws"
// 				if(it->first != "raws" && it->second.idtag != RAWINT_ID)
				if(it->first == bank && it->second.idtag != RAWINT_ID)			
				{	
// 					counter++;
					// integrated raw data
					vector<hitOutput> bankRawHits = getRawIntDataBank(EDT, it->first, &banksMap, verbosity);
					
					// integrated digitized data
					vector<hitOutput> bankDgtHits = getDgtIntDataBank(EDT, it->first, &banksMap, verbosity);
							
					// looping over hits					
					if(bankRawHits.size() == bankDgtHits.size())
					{

						// clean data vector
						layer.clear();
						sector.clear();
						strip.clear();

						pid.clear();
						mpid.clear();
						tid.clear();
						mtid.clear();
						trackE.clear();
						totEdep.clear();
						
						gX.clear();
						gY.clear();
						gZ.clear();
						lX.clear();
						lY.clear();
						lZ.clear();
						
						px.clear();
						py.clear();
						pz.clear();
						vx.clear();
						vy.clear();
						vz.clear();
						mvx.clear();
						mvy.clear();
						mvz.clear();
						t.clear();

						// raw index must be the same as digitized
						for(unsigned h=0; h<bankRawHits.size(); h++)
						{
							// hitn has to be present, has to be same value for raws and dgtz
							if(bankRawHits[h].getIntRawVar("hitn") == bankDgtHits[h].getIntDgtVar("hitn"))
							{
								layer.push_back(bankDgtHits[h].getIntDgtVar("layer"));
								sector.push_back(bankDgtHits[h].getIntDgtVar("sector"));
								strip.push_back(bankDgtHits[h].getIntDgtVar("layer"));
								pid.push_back(bankRawHits[h].getIntRawVar("pid"));
								mpid.push_back(bankRawHits[h].getIntRawVar("mpid"));
								tid.push_back(bankRawHits[h].getIntRawVar("tid"));
								mtid.push_back(bankRawHits[h].getIntRawVar("mtid"));
								trackE.push_back(bankRawHits[h].getIntRawVar("trackE"));
								totEdep.push_back(bankRawHits[h].getIntRawVar("totEdep"));
								gX.push_back(bankRawHits[h].getIntRawVar("<x>"));
								gY.push_back(bankRawHits[h].getIntRawVar("<y>"));
								gZ.push_back(bankRawHits[h].getIntRawVar("<z>"));
								lX.push_back(bankRawHits[h].getIntRawVar("<lx>"));
								lY.push_back(bankRawHits[h].getIntRawVar("<ly>"));
								lZ.push_back(bankRawHits[h].getIntRawVar("<lz>"));
								px.push_back(bankRawHits[h].getIntRawVar("px"));
								py.push_back(bankRawHits[h].getIntRawVar("py"));
								pz.push_back(bankRawHits[h].getIntRawVar("pz"));
								vx.push_back(bankRawHits[h].getIntRawVar("vx"));
								vy.push_back(bankRawHits[h].getIntRawVar("vy"));
								vz.push_back(bankRawHits[h].getIntRawVar("vz"));
								mvx.push_back(bankRawHits[h].getIntRawVar("mvx"));
								mvy.push_back(bankRawHits[h].getIntRawVar("mvy"));
								mvz.push_back(bankRawHits[h].getIntRawVar("mvz"));
								mvz.push_back(bankRawHits[h].getIntRawVar("t"));
 							}
						}
						bankT->Fill();
					}
				}
				
				// generated particle infos is a different tree

				if(it->first == "generated")
				{
					gen_pid.clear();
					gen_p.clear();
					gen_theta.clear();
					gen_phi.clear();
					gen_vx.clear();
					gen_vy.clear();
					gen_vz.clear();
					
					vector<generatedParticle> parts = getGenerated(EDT, getBankFromMap("generated", &banksMap), verbosity);

					for(unsigned p=0; p<parts.size(); p++)
					{
						gen_pid.push_back(parts[p].PID);
						gen_p.push_back(parts[p].momentum.mag());
						gen_theta.push_back(parts[p].momentum.theta()/deg);
						gen_phi.push_back(parts[p].momentum.phi()/deg);
						gen_vx.push_back(parts[p].vertex.x());
						gen_vy.push_back(parts[p].vertex.x());
						gen_vz.push_back(parts[p].vertex.x());
					}
					genT->Fill();
// 					counter++;
				}
// 				cout << counter << endl;			
			}
		}
		
		
		f->Write();
		delete f;
		
	}
	catch (evioException e)
	{
		cerr << e.toString() << endl;
	}
}









