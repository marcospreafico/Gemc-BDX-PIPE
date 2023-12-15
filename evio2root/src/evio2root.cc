// EVIO headers
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"
using namespace evio;

// gemc headers
#include "gemcOptions.h"
#include "gbank.h"
#include "string_utilities.h"


// C++ headers
#include <string>
#include <iostream>
#include <vector>
using namespace std;

// banks header
#include "banks.h"
#include "rootTrees.h"

// disclaimer:
// this is quite bad code, with functions defined in the gemc source code,
// old for loops, and not uniform calls. Some variables are treated a bit differently.
// this is bad but I'm not fixing it 100% because 1: no time and 2: this code will be obsolete when gemc3 will come out

int main(int argc, char **argv)
{
	// argument is input file name.
	try
	{
		// loading options
		goptions gemcOpt;
		gemcOpt.setGoptions();
		
		if(argc == 1)
		{
			char *nargv[2];
			const char buf[] = "-help-all";
			strcpy(nargv[0], buf);
			strcpy(nargv[1], buf);
			gemcOpt.setOptMap(2, nargv);
		}
		
		gemcOpt.setOptMap(argc, argv);
		
		// verbosity
		double verbosity  = gemcOpt.optMap["BANK_VERBOSITY"].arg;
		double MAXN       = gemcOpt.optMap["N"].arg;
		double addEvent   = gemcOpt.optMap["ADDEVN"].arg ;
		double syncRawDgt = gemcOpt.optMap["SYNCRAWDGT"].arg ;

		// list of banks
		string banklist    = gemcOpt.optMap["B"].args ;

// 		vector<string> userBankFileList    = getStringVectorFromStringWithDelimiter(banklist, " ");
// 		vector<string> userBankList    ;
// 		for(auto &bankFile: userBankFileList) {
// 			string bank = bankFile.substr(bankFile.find_last_of("/") + 1);
// 			userBankList.push_back(bank);
// 		}
// 
// 		for(auto &bank: userBankList) {
// 			cout << " Bank found: " << bank << endl;
// 		}

		vector<string> userRawBankList = getStringVectorFromStringWithDelimiter(gemcOpt.optMap["R"].args, " ");
		
		// list of true info variables
		bool writeAll = false;
		vector<string> trueInfoVariables = getStringVectorFromStringWithDelimiter(gemcOpt.optMap["SELECT_RAW_VARIABLES"].args, ":") ;
		if(find(trueInfoVariables.begin(), trueInfoVariables.end(), "all") != trueInfoVariables.end()) {
			writeAll = true;
		}
		
		// for now let's get all systems from TEXT factories, variation "original"
		string factories = "TEXT";
		
		// loading veriables definitions from factories db
		vector<string> whichSystems = get_strings_except(banklist, " all");
		map<string, string> allSystems;
		
		for(unsigned b=0; b<whichSystems.size(); b++) {
			cout << " Bank file to be loaded: " << whichSystems[b] << endl;
			allSystems[whichSystems[b]] = factories;
		}

		map<string, gBank> banksMap = read_banks(gemcOpt, allSystems);

		// selecting input, output file
		string inputfile = gemcOpt.optMap["INPUTF"].args ;
		
		// ROOT file for output; defining tree
		// ROOT does not recognize vector<float> but vector<doubles> are fine
		string outputfile = inputfile.substr(0, inputfile.rfind(".")) + ".root";
		TFile *f = new TFile(outputfile.c_str(),"RECREATE");
		cout << "  > Output File set to: " << outputfile << endl;
		
		
		map<string, rTree> rTrees;
		
		// header bank definitions
		// loading all variables as double
		// maybe we can make it more general later
		rTrees["header"] = rTree(banksMap["header"].bankName, banksMap["header"].bdescription, verbosity);
		for(unsigned i=0; i<banksMap["header"].name.size(); i++) {
			if(banksMap["header"].name[i] == "time")
				rTrees["header"].addVariable(banksMap["header"].name[i], "Ns");
			else
				rTrees["header"].addVariable(banksMap["header"].name[i], "Nd");
		}

		// user bank definitions
		// can't add variables until we see the bank size
		rTrees["userHeader"] = rTree(banksMap["userHeader"].bankName, banksMap["userHeader"].bdescription, verbosity);
		vector<string> userHeaderNames;

		// generated bank definitions
		rTrees["generated"] = rTree(banksMap["generated"].bankName, banksMap["generated"].bdescription, verbosity);
		for(unsigned i=0; i<banksMap["generated"].name.size(); i++) {
			rTrees["generated"].addVariable(banksMap["generated"].name[i], banksMap["generated"].type[i]);
		}
		
		// ancestors bank definitions
		rTrees["ancestors"] = rTree(banksMap["ancestors"].bankName, banksMap["ancestors"].bdescription, verbosity);
		for (unsigned i=0; i<banksMap["ancestors"].name.size(); i++) {
			rTrees["ancestors"].addVariable(banksMap["ancestors"].name[i], banksMap["ancestors"].type[i]);
		}
		
		// hit banks definitions
		for(map<string, gBank>::iterator it = banksMap.begin(); it != banksMap.end(); it++) {

			string thisBankName = it->first;

			bool defineBank = true;

			if(thisBankName == "header")     defineBank = false;
			if(thisBankName == "userHeader") defineBank = false;
			if(thisBankName == "generated")  defineBank = false;
			if(thisBankName == "ancestors")  defineBank = false;
			if(thisBankName == "raws")       defineBank = false;
			if(thisBankName == "psummary")   defineBank = false;

// 			if(find(userBankList.begin(), userBankList.end(), thisBankName) == userBankList.end())   defineBank = false;

			bool defineRawBank = true;
			if(find(userRawBankList.begin(), userRawBankList.end(), thisBankName) == userRawBankList.end()) defineRawBank = false;

			if (defineBank) {

				rTrees[thisBankName] = rTree(banksMap[thisBankName].bankName, banksMap[thisBankName].bdescription, verbosity);

				for(unsigned i=0; i<banksMap[thisBankName].name.size(); i++) {
					// Using Nd for all variables.
					rTrees[thisBankName].addVariable(banksMap[thisBankName].name[i], "Nd");
				}

				// adding raws infos to bank
				if (defineRawBank) {
					for(unsigned i=0; i<banksMap["raws"].name.size(); i++)
						if(banksMap["raws"].name[i] != "hitn") {
							// Using Nd for all variables.
							if(writeAll) {
								rTrees[thisBankName].addVariable(banksMap["raws"].name[i], "Nd");
							} else {
								for(auto varToBeWritten:  trueInfoVariables) {
									if(banksMap["raws"].name[i] == varToBeWritten) {
										rTrees[thisBankName].addVariable(banksMap["raws"].name[i], "Nd");
									}
								}
							}
						}
				}
			}

		}
		
				
		
		// now filling the banks
		// ---------------------
		
		// starting from -3
		// first event is option
		// ++ inside loop
		// equal sign for maxn
		int evn = -3 + addEvent;
		
		evioFileChannel *chan = new evioFileChannel(inputfile, "r", 3000000);
		chan->open();
		
		// skip first event, its the configuration file
		chan->read();
		
		while(chan->read() && (evn++ <= MAXN || MAXN == 0)) {
			evioDOMTree EDT(chan);
			if (evn<-1) continue;			
			
			// read all defined banks
			for(map<string, gBank>::iterator it = banksMap.begin(); it != banksMap.end(); it++) {


				string thisBankName = it->first;


				// header
				if(thisBankName == "header") {
					map<string, double>  headerBank = getHeaderBank(EDT, getBankFromMap("header", &banksMap), 0);
					
					rTrees["header"].init();
					
					for(map<string, double>::iterator head_it = headerBank.begin(); head_it != headerBank.end(); head_it++) {
						// time value is -999
						if(head_it->second == -999) {
							rTrees["header"].insertVariable("time", "Ns", head_it->first);
						} else {
							if(head_it->first == "evn") {
								rTrees["header"].insertVariable(head_it->first, "Nd", head_it->second + addEvent);
							} else {
								rTrees["header"].insertVariable(head_it->first, "Nd", head_it->second);
							}
						}
					}
					rTrees["header"].fill();
				}
				
				// userHeader
				else if (thisBankName == "userHeader") {
					if (userHeaderNames.size() == 0)
					{
						userHeaderNames = getUserHeaderBankNames (EDT, getBankFromMap("userHeader", &banksMap), 0);
						for (vector<string>::iterator it = userHeaderNames.begin(); it != userHeaderNames.end(); it++)
						{
							rTrees["userHeader"].addVariable(*it, "Nd");
						}
					}
					map<string, double>  userHeaderBank = getUserHeaderBank(EDT, getBankFromMap("userHeader", &banksMap), userHeaderNames, 0);
					rTrees["userHeader"].init();
					
					for(map<string, double>::iterator userHead_it = userHeaderBank.begin(); userHead_it != userHeaderBank.end(); userHead_it++) {

						rTrees["userHeader"].insertVariable(userHead_it->first, "Nd", userHead_it->second);
					}
					rTrees["userHeader"].fill();
				}
				
				
				// generated particles
				else if(thisBankName == "generated") {
					vector<generatedParticle> parts = getGenerated(EDT, getBankFromMap(thisBankName, &banksMap), verbosity);
					
					rTrees["generated"].init();
					
					for(unsigned i=0; i<banksMap["generated"].name.size(); i++) {
						string varname = banksMap["generated"].name[i];
						string vartype = banksMap["generated"].type[i];
						
						
						for(unsigned p=0; p<parts.size(); p++) {
							// cout << " p " << p << " " << varname << "  " << vartype << " " << parts[p].getVariableFromStringD(varname) << endl;
							rTrees["generated"].insertVariable(varname, vartype, parts[p].getVariableFromStringI(varname));
							rTrees["generated"].insertVariable(varname, vartype, parts[p].getVariableFromStringD(varname));
						}
					}
					rTrees["generated"].fill();
				}
				
				
				// ancestors
				else if (it->first == "ancestors")
				{
					vector<ancestorInfo> trajs = getAncestors (EDT, getBankFromMap(thisBankName, &banksMap), verbosity);
					
					rTrees["ancestors"].init();
					
					for (unsigned i=0; i<banksMap["ancestors"].name.size(); i++)
					{
						string varname = banksMap["ancestors"].name[i];
						string vartype = banksMap["ancestors"].type[i];


						for(unsigned p=0; p<trajs.size(); p++)
						{
							rTrees["ancestors"].insertVariable(varname, vartype, trajs[p].getVariableFromStringI(varname));
							rTrees["ancestors"].insertVariable(varname, vartype, trajs[p].getVariableFromStringD(varname));
						}
					}
					rTrees["ancestors"].fill();
				}
				
				
				//  hit banks
				else if(thisBankName != "psummary" && thisBankName != "raws") {

					vector<hitOutput> dgtHits = getDgtIntDataBank(EDT, thisBankName, &banksMap, verbosity);
					vector<hitOutput> rawHits = getRawIntDataBank(EDT, thisBankName, &banksMap, verbosity);
					
					// looping over hits
					// hit index must be the same
					unsigned long nrawhits = rawHits.size();
					unsigned long ndgthits = dgtHits.size();

					bool writeRaw = 1;
					if(syncRawDgt == 1 && nrawhits != ndgthits) {
						writeRaw = 0;
					}

					bool fillRawBank = false;

					if(find(userRawBankList.begin(), userRawBankList.end(), thisBankName) != userRawBankList.end()) fillRawBank = true;

					rTrees[thisBankName].init();

					if(writeRaw && fillRawBank) {

						for(unsigned long h=0; h<nrawhits; h++) {
							map<string, double> raws = rawHits[h].getRaws();
							for(map<string, double>::iterator raws_it = raws.begin(); raws_it != raws.end(); raws_it++) {
								if(raws_it->first!= "hitn") {
									rTrees[thisBankName].insertVariable(raws_it->first, "Nd", raws_it->second);
								}
							}
						}
					}

					for(unsigned long h=0; h<ndgthits; h++) {
						map<string, double> dgts = dgtHits[h].getDgtz();
						for(map<string, double>::iterator dgts_it = dgts.begin(); dgts_it != dgts.end(); dgts_it++) {
							rTrees[thisBankName].insertVariable(dgts_it->first, "Nd", dgts_it->second);
						}
					}

					rTrees[thisBankName].filluser();


				}
			}
		}
		chan->close();
		f->Write();
		delete f;
	}
	
	catch (evioException e) {
		cerr << e.toString() << endl;
	}
	
}









