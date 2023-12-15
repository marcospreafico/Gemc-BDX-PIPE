// EVIO headers
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"
using namespace evio;

// gemc headers
#include "options.h"
#include "gbank.h"
#include "string_utilities.h"

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
		
		double verbosity = gemcOpt.optMap["BANK_VERBOSITY"].arg;
		
		
		// list of banks
		string banklist = gemcOpt.optMap["B"].args ;
		
		// for now let's get all systems from TEXT factories, variation "original"
		string factories = "TEXT";
		
		// loading veriables definitions from factories db
		vector<string> whichSystems = get_strings_except(banklist, "all");
		map<string, string> allSystems;
		
		for(unsigned b=0; b<whichSystems.size(); b++)
		allSystems[whichSystems[b]] = factories;
		
		map<string, gBank> banksMap = read_banks(gemcOpt, allSystems);
				
		// selecting input, output file
		string inputfile = gemcOpt.optMap["INPUTF"].args ;
		
		// ROOT file for output; defining tree
		// ROOT does not recognize vector<float> but vector<doubles> are fine
		string outputfile = inputfile.substr(0,inputfile.rfind("."))+".root";
		TFile *f = new TFile(outputfile.c_str(),"RECREATE");
		
		cout << "  > Output File set to: " << outputfile << endl;
		
		// 		int evn;
		vector<int> fluxID;
		vector<double> pid;
		vector<double> mpid;
		vector<double> tid;
		vector<double> mtid;
		vector<double> otid;		
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
		
		// creating the flux tree
		TTree *fluxT = new TTree("fluxT", "flux Tree");
		// 		fluxT->Branch("evn" ,   &evn  );
		
		fluxT->Branch("fluxID", &fluxID);
		
		fluxT->Branch("pid",     &pid);
		fluxT->Branch("mpid",    &mpid);
		fluxT->Branch("tid",     &tid);
		fluxT->Branch("mtid",    &mtid);
		fluxT->Branch("otid",    &otid);
		
		fluxT->Branch("trackE",  &trackE);
		fluxT->Branch("totEdep", &totEdep);
		
		fluxT->Branch("gX" ,     &gX);
		fluxT->Branch("gY" ,     &gY);
		fluxT->Branch("gZ" ,     &gZ);
		fluxT->Branch("lX" ,     &lX);
		fluxT->Branch("lY" ,     &lY);
		fluxT->Branch("lZ" ,     &lZ);
		
		fluxT->Branch("px" ,     &px);
		fluxT->Branch("py" ,     &py);
		fluxT->Branch("pz" ,     &pz);
		fluxT->Branch("vx" ,     &vx);
		fluxT->Branch("vy" ,     &vy);
		fluxT->Branch("vz" ,     &vz);
		fluxT->Branch("mvx" ,    &mvx);
		fluxT->Branch("mvy" ,    &mvy);
		fluxT->Branch("mvz" ,    &mvz);
		fluxT->Branch("t" ,      &t);
		
		// creating the generated particle info tree
		vector<double> gen_pid;
		vector<double> gen_p;
		vector<double> gen_theta;
		vector<double> gen_phi;
		vector<double> gen_px;
		vector<double> gen_py;
		vector<double> gen_pz;		
		vector<double> gen_vx;
		vector<double> gen_vy;
		vector<double> gen_vz;
		vector<string> gen_dname;
		vector<int> gen_stat;
		vector<double> gen_etot;
		vector<double> gen_time;		
		
		
		TTree *genT = new TTree("genT", "generated particles Tree");
		
		genT->Branch("pid"  ,   &gen_pid);
		genT->Branch("p"    ,   &gen_p);
		genT->Branch("theta",   &gen_theta);
		genT->Branch("phi"  ,   &gen_phi);
		genT->Branch("px"   ,   &gen_px);
		genT->Branch("py"   ,   &gen_py);
		genT->Branch("pz"   ,   &gen_pz);		
		genT->Branch("vx"   ,   &gen_vx);
		genT->Branch("vy"   ,   &gen_vy);
		genT->Branch("vz"   ,   &gen_vz);
		genT->Branch("dname"   ,   &gen_dname);
		genT->Branch("stat"   ,   &gen_stat);
		genT->Branch("etot"   ,   &gen_etot);
		genT->Branch("time"   ,   &gen_time);		
		
		// creating the header info tree
		int evn,ngen,evn_typ;
		double var1,var2,var3,var4,var5,var6,var7,var8,var9;
		
		TTree *headerT = new TTree("headerT", "header Tree");
		headerT->Branch("evn" ,   &evn);
		
		headerT->Branch("ngen" ,   &ngen);
		headerT->Branch("evn_typ" ,   &evn_typ);
		
		headerT->Branch("var1"  ,   &var1);
		headerT->Branch("var2"  ,   &var2);
		headerT->Branch("var3"  ,   &var3);
		headerT->Branch("var4"  ,   &var4);
		headerT->Branch("var5"  ,   &var5);
		headerT->Branch("var6"  ,   &var6);
		headerT->Branch("var7"  ,   &var7);
		headerT->Branch("var8"  ,   &var8);
		headerT->Branch("var9"  ,   &var9);
		
		evioFileChannel *chan = new evioFileChannel(inputfile, "r", 3000000);
		chan->open();
		
		int counter_flux=0,counter_gen=0,counter_header=0;
		while(chan->read())
		{
			evioDOMTree EDT(chan);
			
			// read hit banks
			for(map<string, gBank>::iterator it = banksMap.begin(); it != banksMap.end(); it++)
			{
				counter_flux++;
				if (counter_flux==1) continue;
				// only considering data banks (bankType.size() > 0) different from "raws"
				// 				if(it->first != "raws" && it->second.idtag != RAWINT_ID)
				if(it->first == "flux" && it->second.idtag != RAWINT_ID)
				{
					// integrated raw data
					vector<hitOutput> fluxRawHits = getRawIntDataBank(EDT, it->first, &banksMap, verbosity);
					
					// integrated digitized data
					vector<hitOutput> fluxDgtHits = getDgtIntDataBank(EDT, it->first, &banksMap, verbosity);
					
					// looping over hits
					if(fluxRawHits.size() == fluxDgtHits.size())
					{
						
						// clean data vector
						fluxID.clear();
						
						pid.clear();
						mpid.clear();
						tid.clear();					
						mtid.clear();
						otid.clear();
						
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
						for(unsigned h=0; h<fluxRawHits.size(); h++)
						{
							// hitn has to be present, has to be same value for raws and dgtz
							if(fluxRawHits[h].getIntRawVar("hitn") == fluxDgtHits[h].getIntDgtVar("hitn"))
							{
								fluxID.push_back(fluxDgtHits[h].getIntDgtVar("id"));
								
								pid.push_back(fluxRawHits[h].getIntRawVar("pid"));
								mpid.push_back(fluxRawHits[h].getIntRawVar("mpid"));
								tid.push_back(fluxRawHits[h].getIntRawVar("tid"));
								mtid.push_back(fluxRawHits[h].getIntRawVar("mtid"));
								otid.push_back(fluxRawHits[h].getIntRawVar("otid"));
							
								trackE.push_back(fluxRawHits[h].getIntRawVar("trackE"));							totEdep.push_back(fluxRawHits[h].getIntRawVar("totEdep"));
								gX.push_back(fluxRawHits[h].getIntRawVar("<x>"));
								gY.push_back(fluxRawHits[h].getIntRawVar("<y>"));
								gZ.push_back(fluxRawHits[h].getIntRawVar("<z>"));
								lX.push_back(fluxRawHits[h].getIntRawVar("<lx>"));
								lY.push_back(fluxRawHits[h].getIntRawVar("<ly>"));
								lZ.push_back(fluxRawHits[h].getIntRawVar("<lz>"));
								px.push_back(fluxRawHits[h].getIntRawVar("px"));
								py.push_back(fluxRawHits[h].getIntRawVar("py"));
								pz.push_back(fluxRawHits[h].getIntRawVar("pz"));
								vx.push_back(fluxRawHits[h].getIntRawVar("vx"));
								vy.push_back(fluxRawHits[h].getIntRawVar("vy"));
								vz.push_back(fluxRawHits[h].getIntRawVar("vz"));
								mvx.push_back(fluxRawHits[h].getIntRawVar("mvx"));
								mvy.push_back(fluxRawHits[h].getIntRawVar("mvy"));
								mvz.push_back(fluxRawHits[h].getIntRawVar("mvz"));
								t.push_back(fluxRawHits[h].getIntRawVar("<t>"));
								
//   cout << fluxDgtHits[h].getIntDgtVar("id") << " " << fluxRawHits[h].getIntRawVar("pid")  << " " << fluxRawHits[h].getIntRawVar("mpid")  << " " << fluxRawHits[h].getIntRawVar("tid")  << " " << fluxRawHits[h].getIntRawVar("mtid")  << " " << fluxRawHits[h].getIntRawVar("otid")  << " " << fluxRawHits[h].getIntRawVar("trackE")  << " " << fluxRawHits[h].getIntRawVar("totEdep")  << " " << fluxRawHits[h].getIntRawVar("<x>")  << " " << fluxRawHits[h].getIntRawVar("<y>")  << " " << fluxRawHits[h].getIntRawVar("<z>")  << " " << fluxRawHits[h].getIntRawVar("<lx>")  << " " << fluxRawHits[h].getIntRawVar("<ly>")  << " " << fluxRawHits[h].getIntRawVar("<lz>")  << " " << fluxRawHits[h].getIntRawVar("px")   << " " << fluxRawHits[h].getIntRawVar("py")   << " " << fluxRawHits[h].getIntRawVar("pz")   << " " << fluxRawHits[h].getIntRawVar("vx")   << " " << fluxRawHits[h].getIntRawVar("vy")   << " " << fluxRawHits[h].getIntRawVar("vz")   << " " << fluxRawHits[h].getIntRawVar("mvx")   << " " << fluxRawHits[h].getIntRawVar("mvy")   << " " << fluxRawHits[h].getIntRawVar("mvz")   << " " << fluxRawHits[h].getIntRawVar("<t>")   << endl;
 							}
						}
						fluxT->Fill();
					}
				}
				
				// generated particle infos is a different tree
				
				if(it->first == "generated")
				{
					counter_gen++;
					if (counter_gen==1) continue;
					
					gen_pid.clear();
					gen_p.clear();
					gen_theta.clear();
					gen_phi.clear();
					gen_px.clear();
					gen_py.clear();
					gen_pz.clear();					
					gen_vx.clear();
					gen_vy.clear();
					gen_vz.clear();
					gen_dname.clear();
					gen_stat.clear();
					gen_etot.clear();
					gen_time.clear();					
					
					vector<generatedParticle> parts = getGenerated(EDT, getBankFromMap("generated", &banksMap), verbosity);
					
					for(unsigned p=0; p<parts.size(); p++)
					{
						gen_pid.push_back(parts[p].PID);
						gen_p.push_back(parts[p].momentum.mag()/MeV);
						gen_theta.push_back(parts[p].momentum.theta()/deg);
						gen_phi.push_back(parts[p].momentum.phi()/deg);
						gen_px.push_back(parts[p].momentum.x()/MeV);
						gen_py.push_back(parts[p].momentum.y()/MeV);
						gen_pz.push_back(parts[p].momentum.z()/MeV);
						gen_vx.push_back(parts[p].vertex.x()*10.);
						gen_vy.push_back(parts[p].vertex.y()*10.);
						gen_vz.push_back(parts[p].vertex.z()*10.);
// 						gen_dname.push_back(parts[p].dname);
// 						gen_stat.push_back(parts[p].stat);
// 						gen_etot.push_back(parts[p].etot/MeV);
// 						gen_time.push_back(parts[p].time);					
// 						cout << parts[p].PID << " " << parts[p].momentum.mag()/MeV << " " << parts[p].momentum.theta()/deg << " " << parts[p].momentum.phi()/deg << " " << parts[p].momentum.x()/MeV << " " << parts[p].momentum.y()/MeV << " " << parts[p].momentum.z()/MeV << " " << parts[p].vertex.x()*10. << " " << parts[p].vertex.y()*10. << " " << parts[p].vertex.z()*10. << endl;
					}
					genT->Fill();
				}
				
				if(it->first == "header")
				{
					counter_header++;
					if (counter_header==1) continue;
					
					map<string, double>  headerBank = getHeaderBank(EDT, getBankFromMap("header", &banksMap), verbosity);
					
					for(map<string, double> :: iterator itit = headerBank.begin(); itit != headerBank.end(); itit++)
					{
						// 	cout << getBankFromMap("header", &banksMap).getVarId(itit->first) << " " << itit->first << ":\t" << itit->second << endl;
						switch (getBankFromMap("header", &banksMap).getVarId(itit->first)) {
								case 2:
								evn=itit->second;
								break;
								case 3:
								evn_typ=itit->second;
								break;
								case 4:
								var4=itit->second;
								break;
								case 5:
								var3=itit->second;
								break;
								case 6:
								ngen=itit->second;
								break;
								case 7:
								var1=itit->second;
								break;
								case 8:
								var2=itit->second;
								break;
								case 9:
								var5=itit->second;
								break;
								case 10:
								var6=itit->second;
								break;
								case 11:
								var7=itit->second;
								break;
								case 12:
								var8=itit->second;
								break;
								case 13:
								var9=itit->second;
								break;
							default:
								break;
						}
					}
					headerT->Fill();
					// 	cout << evn << "\t" << evn_typ << "\t" << ngen << "\t" << var1 << "\t" << var2 << "\t" << var3 << "\t" << var4 << "\t" << var5 << "\t" << var6 << "\t" << var7 << "\t" << var8 << "\t" << var9 << endl;
					
				}
				
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









