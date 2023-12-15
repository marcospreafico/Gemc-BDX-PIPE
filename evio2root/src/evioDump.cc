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
#include <set>
using namespace std;

// banks header
#include "banks.h"


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
		
		// reading input
		string inputfile = gemcOpt.optMap["INPUTF"].args ;
		evioFileChannel *chan = new evioFileChannel(inputfile, "r", 3000000);
		chan->open();
		while(chan->read())
		{
			evioDOMTree EDT(chan);
			
			// get simulation conditions
			readConditions(EDT, verbosity);
	
			// get header bank
			map<string, double> headerBank = getHeaderBank(EDT, getBankFromMap("header", &banksMap), verbosity);

			cout << endl;
		
			for(map<string, gBank>::iterator it = banksMap.begin(); it != banksMap.end(); it++)
			{
				// integrated raw data
				getRawIntDataBank(EDT, it->first, &banksMap, verbosity);
				
				// integrated digitized data
				getDgtIntDataBank(EDT, it->first, &banksMap, verbosity);
			}
		}
		
		
	}
	catch (evioException e)
	{
		cerr << e.toString() << endl;
	}
}
